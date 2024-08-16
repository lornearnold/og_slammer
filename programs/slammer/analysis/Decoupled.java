/* This file is in the public domain. */
// 1 - PREVIOUS ITERATION
// 2 - CURRENT ITERATION

package slammer.analysis;

import slammer.*;

public class Decoupled extends DeCoupledCommon
{
	private double deltacc;
	private double sdot[];

	private double time;

	public double Decoupled(double[] ain_p, double uwgt_p, double height_p, double vs_p, double damp1_p, double refstrain_p, double dt_p, double scal_p, double g_p, double vr_p, double[][] ca, boolean dv3_p)
	{
		// assign all passed parameters to the local data
		uwgt = uwgt_p; // UNIT WEIGHT? (MAYBE SPECIFIC WEIGHT)
		height = height_p; // HEIGHT OF SLOPE
		vs = vs_p; // SHEAR WAVE VELOCITY ABOVE SURFACE (IN SLIDING MASS)
		vs1 = vs; // COPY OF vs INTO 'PREVIOUS' ITERATION
		damp1 = damp1_p; // DAMPING RATIO
		damp = damp1; // COPY OF damp1
		dt = dt_p; // DT
		scal = scal_p; // SCALE FACTOR?
		g = g_p; 
		vr = vr_p; // SHEAR WAVE VELOCITY BELOW SURFACE (IN BASE)
		dv3 = dv3_p; // EQUIVALENT LINEAR (IF TRUE) OR LINEAR ELASTIC (IF FALSE)
		ain = ain_p; // GROUND MOTION ACCELERATION VALUES

		if ((g < 33.0) && (g > 32.0))
			uwgt = 120.0;
		else
			uwgt = 20.0;

		// init graphing
		setValueSize(dt);

		// copy ca into disp and mu
		// mu AND disp RELATED TO CRITICAL ACCELERATION. LENGTH > 1 IF CRITICAL ACCELERATION VARIES WITH DISPLACEMENT
		disp = new double[ca.length];
		mu = new double[ca.length];
		for(int i = 0; i < ca.length; i++)
		{
			disp[i] = ca[i][0];
			mu[i] = ca[i][1];
		}

		nmu = ca.length; // LENGTH OF mu TABLE

		npts = ain.length; // NUMBER OF POINTS IN GROUND MOTION

		avgacc = new double[npts];
		s = new double[npts];
		sdot = new double[npts];
		u = new double[npts];
		udot = new double[npts];
		udotdot = new double[npts];

		deltacc = 0.0;
		// u1, udot, udotdot1 GROUND DISP, VELOCITY, ACCELERATION RESPECTIVELY
		u1 = 0.0;
		udot1 = 0.0;
		udotdot1 = 0.0;

		acc1 = 0.0;
		acc2 = 0.0;
		mx = 0.0;
		mx1 = 0.0;
		mmax = 0.0;
		gameff1 = 0.0;
		n = 100.0;
		o = 100.0;

		rho = uwgt / g; // DENSITY

		dampf = 55.016 * Math.pow((vr / vs), -0.9904) / 100.0;
		if(dampf > 0.2)
			dampf = 0.2;

		scal *= -1.0;

		// for each mode calculate constants for Slammer algorithm

		beta = 0.25;
		gamma = 0.5;
		Mtot = rho * height;
		slide = true;
		qq = 1;

		omega = Math.PI * vs / (2.0 * height);
		// L AND M CONSTANTS FROM RATHJE 1999?
		L = -2.0 * rho * height / Math.PI * Math.cos(Math.PI);
		M = rho * height / 2.0;

		n = 100.0;
		o = 100.0;
		gamref = refstrain_p;
		damp = damp1 + dampf;

		// loop for time steps in time histories

		// for equivalent linear
		if(dv3)
			eq();

		omega = Math.PI * vs / (2.0 * height);

		// calculate final dynamic response using original prop's for LE analysis and EQL properties for EQL analysis
		for(j = 1; j <= npts; j++)
		{
			d_setupstate();
			d_response();
		}

		slide = false;
		time = 0.0;

		avg_acc();

		_kmax = mmax;
		_vs = vs;
		_damp = damp;
		_dampf = dampf;
		_omega = omega;

		// calculate decoupled displacements
		for(j = 1; j <= npts; j++)
		{
			d_sliding();
			store(Math.abs(s[j - 1]));
			residual_mu();
		}

		end(Math.abs(s[npts - 1]));
		return Math.abs(s[npts - 1]);
	}

	private void d_sliding()
	{
		// calculate decoupled displacements

		double deltacc;

		if(j == 1)
			deltacc = avgacc[j - 1];
		else
			deltacc = avgacc[j - 1] - avgacc[j - 2];

		if(j == 1) // added
		{
			sdot[j - 1] = 0;
			s[j - 1] = 0;
		}
		else if(!slide)
		{
			sdot[j - 1] = 0;
			s[j - 1] = s[j - 2];
		}
		else
		{
			sdot[j - 1] = sdot[j - 2] + (mu[qq - 1] * g - avgacc[j - 2]) * dt - 0.5 * deltacc * dt;
			s[j - 1] = s[j - 2] - sdot[j - 2] * dt - 0.5 * dt * dt * (mu[qq - 1] * g - avgacc[j - 2]) + deltacc * dt * dt / 6.0;
		}

		if(!slide)
		{
			if(avgacc[j - 1] > mu[qq - 1] * g)
				slide = true;
		}
		else
		{
			if(sdot[j - 1] >= 0.0)
			{
				slide = false;

				if(j == 1)
					s[j - 1] = 0;
				else
					s[j - 1] = s[j - 2];

				sdot[j - 1] = 0.0;
			}
		}
	}
}
