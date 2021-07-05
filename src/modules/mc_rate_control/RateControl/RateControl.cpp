/****************************************************************************
 *
 *   Copyright (c) 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file RateControl.cpp
 */

#include <RateControl.hpp>
#include <px4_platform_common/defines.h>

/** Begin SYSU Lab

 * Used for Model Based Control.

 * Model Based Controller outputs commands of thrust and torques around three-axis,

 * while the commands send by px4 to low-level mixer and motors are normalized values in 0~1(thrust) or -1~1(torque).

 * By converting F&T commands to normalized commands, we can implement Model Based Control methods readily,

 * using existing mixer and saturation limiter.

 */

/* Model of QAV250 */

#define HALF_LENGTH	0.101f

#define HALF_WIDTH	0.080f


#define C_M		0.0097f		// C_M = drag torque / motor thrust.


#define I_XX		3e-3f

#define I_YY		2.7e-3f

#define I_ZZ		6.0e-3f

#define GRA_ACC		9.8066f

#define MAV_MASS	0.703f
// End SYSU Lab


using namespace matrix;

void RateControl::setGains(const Vector3f &P, const Vector3f &I, const Vector3f &D)
{
	_gain_p = P;
	_gain_i = I;
	_gain_d = D;
}

void RateControl::setSaturationStatus(const MultirotorMixer::saturation_status &status)
{
	_mixer_saturation_positive[0] = status.flags.roll_pos;
	_mixer_saturation_positive[1] = status.flags.pitch_pos;
	_mixer_saturation_positive[2] = status.flags.yaw_pos;
	_mixer_saturation_negative[0] = status.flags.roll_neg;
	_mixer_saturation_negative[1] = status.flags.pitch_neg;
	_mixer_saturation_negative[2] = status.flags.yaw_neg;
}

Vector3f RateControl::update(const Vector3f &rate, const Vector3f &rate_sp, const Vector3f &angular_accel,
			     const float dt, const bool landed)
{
	// angular rates error
	Vector3f rate_error = rate_sp - rate;

	// PID control with feed forward
	//const Vector3f torque = _gain_p.emult(rate_error) + _rate_int - _gain_d.emult(angular_accel) + _gain_ff.emult(rate_sp);

	//Begin SYSU Lab
	/* using PD controller for linear part */

	Vector3f att_ctrl = _gain_p.emult(rate_error) - _gain_d.emult(angular_accel);
	/* Compensate nonlinear gyro moment and ... */

	Vector3f torque_affix, compen_torque, compen_control;
	torque_affix(0) = (I_ZZ - I_YY) * rate(1) * rate(2);
	torque_affix(1) = (I_XX - I_ZZ) * rate(0) * rate(2);
	torque_affix(2) = (I_YY - I_XX) * rate(0) * rate(1);
	compen_torque = torque_affix;

	/* Convert torque command to normalized input */
	compen_control = torque_to_attctrl(compen_torque);

	/* Plus all control as output*/
	const Vector3f torque = att_ctrl + compen_control;

	//End SYSU Lab

	// update integral only if we are not landed
	if (!landed) {
		updateIntegral(rate_error, dt);
	}

	return torque;
}

void RateControl::updateIntegral(Vector3f &rate_error, const float dt)
{
	for (int i = 0; i < 3; i++) {
		// prevent further positive control saturation
		if (_mixer_saturation_positive[i]) {
			rate_error(i) = math::min(rate_error(i), 0.f);
		}

		// prevent further negative control saturation
		if (_mixer_saturation_negative[i]) {
			rate_error(i) = math::max(rate_error(i), 0.f);
		}

		// I term factor: reduce the I gain with increasing rate error.
		// This counteracts a non-linear effect where the integral builds up quickly upon a large setpoint
		// change (noticeable in a bounce-back effect after a flip).
		// The formula leads to a gradual decrease w/o steps, while only affecting the cases where it should:
		// with the parameter set to 400 degrees, up to 100 deg rate error, i_factor is almost 1 (having no effect),
		// and up to 200 deg error leads to <25% reduction of I.
		float i_factor = rate_error(i) / math::radians(400.f);
		i_factor = math::max(0.0f, 1.f - i_factor * i_factor);

		// Perform the integration using a first order method
		float rate_i = _rate_int(i) + i_factor * _gain_i(i) * rate_error(i) * dt;

		// do not propagate the result if out of range or invalid
		if (PX4_ISFINITE(rate_i)) {
			_rate_int(i) = math::constrain(rate_i, -_lim_int(i), _lim_int(i));
		}
	}
}

void RateControl::getRateControlStatus(rate_ctrl_status_s &rate_ctrl_status)
{
	rate_ctrl_status.rollspeed_integ = _rate_int(0);
	rate_ctrl_status.pitchspeed_integ = _rate_int(1);
	rate_ctrl_status.yawspeed_integ = _rate_int(2);
}

/* Added by SYSU Lab: Zhiwei Hou
* Transform computed torques to equivalent normalized inputs
*/

Vector3f RateControl::torque_to_attctrl(matrix::Vector3f &computed_torque)

{

	/* motor maximum thrust model */
	//float thrust_max = 5.488f * sinf(_battery_status.voltage_filtered_v * 0.4502f + 2.2241f);
	float thrust_max = 5.488f;
	/* mixer matrix */
	float array_mixer[4][3] = {
			{-0.707107f,  0.707107f,  1.0f},
			{0.707107f,  -0.707107f,  1.0f},
			{0.707107f,  0.707107f,  -1.0f},
			{-0.707107f,  -0.707107f,  -1.0f}
	};
	Matrix<float, 4, 3> mixer_matrix(array_mixer);
	/* matrix from thrust of four propellers to torque about 3-axes */
	float array_torque[3][4] = {
			{-HALF_LENGTH, HALF_LENGTH, HALF_LENGTH, -HALF_LENGTH},
			{HALF_WIDTH, -HALF_WIDTH, HALF_WIDTH, -HALF_WIDTH},
			{C_M, C_M, -C_M, -C_M}
	};
	Matrix<float, 3, 4> gentrq_matrix(array_torque);
	/* input = (Gamma * Mixer * Tmax)^(-1) * computed_torque */

	SquareMatrix<float, 3> gentrq_mixer = Matrix<float, 3, 3>(gentrq_matrix * mixer_matrix);
	Matrix<float, 3, 3> trq_to_attctrl = gentrq_mixer.I() / thrust_max;
	//PX4_INFO("trq_to_attctrl's diag = %d, %d, %d", (int)(trq_to_attctrl(0,0)*1000.0f), (int)(trq_to_attctrl(1,1)*1000.0f), (int)(trq_to_attctrl(2,2)*1000.0f));
	return trq_to_attctrl * computed_torque;
}




