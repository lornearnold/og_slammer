/*
 * Analysis.java - the scientific analysis algorithms
 *
 * Copyright (C) 2002 Matthew Jibson (dolmant@dolmant.net)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

/* $Id: Analysis.java,v 1.4 2003/08/07 02:41:54 dolmant Exp $ */

package newmark;

import java.text.DecimalFormat;
import java.io.*;
import org.jfree.data.xy.XYSeries;
import newmark.gui.*;

public class Analysis
{
	public static final String fmtFour              = "0.0000";
	public static final String fmtThree             = "0.000";
	public static final String fmtTwo               = "0.00";
	public static final String fmtOne               = "0.0";
	public static final String fmtZero              = "0";

	public static final double Gcmss	= 980.665;
	public static final double Gftss  = Gcmss * 0.032808399;

	public static final int each = 1;
	private static double time;
	private static int eachAt;
	private static double dint;
	public static XYSeries xys;
	private static double last;
	private static boolean skipped;

	private static int perSec = 5;
	private static double interval = 1.0 / (double)perSec;
	private static double timeStor;

	private static void setValueSize(final double Dint)
	{
		time = 0;
		eachAt = each;
		dint = Dint;
		xys = new XYSeries("");
		last = -1;
		skipped = false;
		timeStor = 0;
	}

	private static void store(final double d)
	{
		if(d == last)
		{
			skipped = true;
		}
		else
		{
			if(skipped)
			{
				realStore(last, time - dint);
				skipped = false;
			}

			if(time >= (timeStor + interval))
			{
				realStore(d, time);
				timeStor = time;
			}
		}

		time += dint;
	}

	private static void end(final double d)
	{
		if(skipped)
			realStore(last, time - dint);
		realStore(d, time);
	}

	private static void realStore(final double d, final double time)
	{
		try {xys.add(new Float(time), new Float(d));}
		catch (Exception e) {}
		last = d;
	}

	/* standard functions */

	public static double log10(final double val)
	{
		return Math.log(val) / Math.log(10.0);
	}

	public static double sign(final double val)
	{
		if(val >= 0)
			return 1;
		else
			return -1;
	}

	/* the algorithms/programs
	 *
	 * Note: the arias, dobry, redigitize, and rigorous newmark algorithms were
	 * originally written by Ray Wilson in GW Basic, and then ported to C++ and
	 * Java.
	 */

	public static String Arias(DoubleList data, final double di)
	{
		DecimalFormat fmt = new DecimalFormat(fmtThree);
		return fmt.format(AriasDobry(data, di, false));
	}

	/* Arias, with optional dobry.  If the boolean dobry is true, this function
	 * will return the dobry duration as opposed to the arias intensity.
	 */
	private static double AriasDobry(DoubleList data, final double di, boolean dobry)
	{
		Double val;
		double d, u = 0, a, x, j, z, b, t, k = 0, res;
		d = di;
		j = .1603 * d * u;
		data.reset();
		while((val = data.each()) != null)
		{
			a = val.doubleValue();
			a /= 100;
			x = a * a;
			u += x;
		}
		j = .1603 * d * u;
		if(dobry)
		{	z = .95 * j / (.1603 * d);
			b = .05 * j / (.1603 * d);
			u = 0;
			data.reset();
			while((val = data.each()) != null)
			{
				a = val.doubleValue();
				a /= 100;
				x = a * a;
				u += x;
				if(u > z || u < b) k--;
				k++;
			}
			t = k * d;
			return t;
		}
		else
		{
			return j;
		}
	}

	public static String CM_GS(DoubleList data, FileWriter ofile) throws IOException
	{
		final double val = 1.0 / Gcmss;
		return Mult(data, ofile, val);
	}

	public static String Count(DoubleList data)
	{
		DecimalFormat fmt = new DecimalFormat(fmtZero);
		return fmt.format(data.size());
	}

	public static String Dobry(DoubleList data, final double di)
	{
		DecimalFormat fmt = new DecimalFormat(fmtOne);
		return fmt.format(AriasDobry(data, di, true));
	}

	public static String PGA(DoubleList data)
	{
		DecimalFormat fmt = new DecimalFormat(fmtThree);
		return fmt.format(FindMax(data) / Gcmss); // store in g's, but expect to be in cm/s/s
	}

	private static double FindMax(DoubleList data)
	{
		Double val;
		double here, max = 0;
		data.reset();
		while((val = data.each()) != null)
		{
			here = Math.abs(val.doubleValue());
			if(here < 0) here = -here;
			if(here > max) max = here;
		}
		return max;
	}

	public static String GS_CM(DoubleList data, FileWriter ofile) throws IOException
	{
		final double val = Gcmss;
		return Mult(data, ofile, val);
	}

	public static String Mult(DoubleList data, FileWriter ofile, final double value) throws IOException
	{
		Double val;
		double temp;
		data.reset();
		while((val = data.each()) != null)
		{
			temp = val.doubleValue();
			temp *= value;
			ofile.write(Double.toString(temp));
			ofile.write('\n');
		}
		ofile.close();
		return null;
	}

	public static String Peapick(DoubleList data, FileWriter ofile, final double value) throws IOException
	{
		Double val;
		int top = data.size() - 1;
		int indexl = 0, index = 0, indexr = top;

		data.reset();
		while((val = data.each()) != null)
		{
			if(Math.abs(val.doubleValue()) >= value)
			{
				indexl = index - 50;
				if(indexl < 0)
					indexl = 0;
				break;
			}
			index++;
		}

		data.end();
		index = top;
		while((val = data.eachP()) != null)
		{
			if(Math.abs(val.doubleValue()) >= value)
			{
				indexr = index + 50;
				if(indexr > top)
					indexr = top;
				break;
			}
			index--;
		}

		// put data.current at indexl
		data.reset();
		for(int i = 0; i < indexl; i++)
			data.next();

		// start writing to the file
		for(int i = indexl; i <= indexr; i++)
		{
			ofile.write(data.each().toString());
			ofile.write('\n');
		}

		ofile.close();

		return null;
	}

	public static String NewmarkRigorousDual(DoubleList data, final double d, final double ca, final double ta, final double mult)
	{
		Double val;
		double a, n, q=0, r=0, s=0, y=0, v=0, u=0;

		final double l = Math.toRadians(ta);
		final double g = Math.sin(l) * Gcmss;
		final double t = (ca * Gcmss) + g;

		setValueSize(d);
		data.reset();
		while((val = data.each()) != null)
		{
			a = (val.doubleValue() * mult) + g;

			if(a == 0)
			{
				store(u);
				continue;
			}

			if(Math.abs(v) < .0001)
			{
				if(Math.abs(a) > t)
				{
					n = sign(a);
				}
				else
				{
					n = a / t;
				}
			}
			else
			{
				n = sign(v);
			}
			y = a - n * t;
			v = r + d / 2.0 * (y + s);
			if (!(r == 0.0 || (v / r) > 0))
			{
				v = 0;
				y = 0;
			}
			u = q + d / 2.0 * (v + r);
			q = u;
			r = v;
			s = y;
			if(mult > 0.0) store(u);
		}

		DecimalFormat fmt = new DecimalFormat(fmtOne);
		end(u);
		return fmt.format(u);
	}

	public static String NewmarkRigorous(DoubleList data, final double di, final double ca, final double mult)
	{
		Double val;
		double t, d, q = 0, r = 0, s = 0, y = 0, v = 0, u = 0, a, n;
		t = ca;
		t *= Gcmss;
		d = di;
		setValueSize(di);
		data.reset();
		int count = 0;
		while((val = data.each()) != null)
		{
			a = val.doubleValue() * mult;
			if(a == 0)
			{
				store(u);
				continue;
			}
			if(v < .0001)
			{
				if(Math.abs(a) > t )
				{
					n = sign(a);
				}
				else
				{
					n = a / t;
				}
			}
			else
			{
				n = 1;
			}
			y = a - n * t;
			v = r + d / 2.0 * (y + s);
			if (v <= 0.0)
			{
				v = 0;
				y = 0;
			}
			u = q + d / 2.0 * (v + r);
			q = u;
			r = v;
			s = y;
			if(mult > 0.0) store(u);
		}

		DecimalFormat fmt = new DecimalFormat(fmtOne);
		end(u);
		return fmt.format(u);
	}

	public static String NewmarkRigorousDisp(DoubleList data, final double di, final double[][] disp, final double mult)
	{
		Double val;
		double t, d, q = 0, r = 0, s = 0, y = 0, v = 0, u = 0, a, n;
		t = disp[0][1];
		t *= Gcmss;
		int pos = 0;
		double prop;
		d = di;
		setValueSize(di);
		data.reset();
		while((val = data.each()) != null)
		{
			a = val.doubleValue() * mult;
			if(a == 0)
			{
				store(u);
				continue;
			}
			if(v < .0001)
			{
				if(Math.abs(a) > t )
					n = sign(a);
				else
					n = a / t;
			}
			else
				n = 1;
			y = a - n * t;
			v = r + d / 2 * (y + s);
			if (v <= 0)
			{
				v = 0;
				y = 0;
			}
			u = q + d / 2 * (v + r);
			q = u;
			r = v;
			s = y;
			if(mult > 0) store(u);
			if(pos == disp.length - 1)
				continue;
			while(u > disp[pos + 1][0])
			{
				pos++;
				if(pos == disp.length - 1) break;
			}
			if(pos == disp.length - 1)
			{
				t = Gcmss * disp[pos][1];
				continue;
			}
			prop = (u - disp[pos][0]) / (disp[pos + 1][0] - disp[pos][0]);
			t = Gcmss * (disp[pos][1] - (disp[pos][1] - disp[pos + 1][1]) * prop);
		}

		DecimalFormat fmt = new DecimalFormat(fmtOne);
		end(u);
		return fmt.format(u);
	}

	public static String NewmarkRigorousTime(DoubleList data, final double di, final double[][] disp, final double mult)
	{
		Double val;
		double t, d, q = 0, r = 0, s = 0, y = 0, v = 0, u = 0, a, n;
		double time = 0;
		t = disp[0][1];
		t *= Gcmss;
		int pos = 0;
		double prop;
		setValueSize(di);
		d = di;
		data.reset();
		while((val = data.each()) != null)
		{
			a = val.doubleValue() * mult;
			if(a == 0)
			{
				store(u);
				continue;
			}
			if(v < .0001)
			{
				if(Math.abs(a) > t )
					n = sign(a);
				else
					n = a / t;
			}
			else
				n = 1;
			y = a - n * t;
			v = r + d / 2 * (y + s);
			if (v <= 0)
			{
				v = 0;
				y = 0;
			}
			u = q + d / 2 * (v + r);
			q = u;
			r = v;
			s = y;
			if(mult > 0) store(u);
			if(pos == disp.length - 1)
				continue;
			while(time > disp[pos + 1][0])
			{
				pos++;
				if(pos == disp.length - 1) break;
			}
			if(pos == disp.length - 1)
			{
				t = Gcmss * disp[pos][1];
				continue;
			}
			prop = (time - disp[pos][0]) / (disp[pos + 1][0] - disp[pos][0]);
			t = Gcmss * (disp[pos][1] - (disp[pos][1] - disp[pos + 1][1]) * prop);
			time += d;
		}

		DecimalFormat fmt = new DecimalFormat(fmtOne);
		end(u);
		return fmt.format(u);
	}

	public static String JibsonAndOthers(final double arias, final double ca)
	{
		DecimalFormat fmt = new DecimalFormat(fmtOne);
		return fmt.format(Math.pow(10, 1.521 * log10(arias) - 1.993 * log10(ca) -1.546));
	}

	public static String AmbraseysAndMenu(final double pga, final double ca)
	{
		final double ratio = ca / pga;
		DecimalFormat fmt = new DecimalFormat(fmtOne);
		return fmt.format(Math.pow(10, 0.90 + log10(Math.pow(1.0 - ratio, 2.53) * Math.pow(ratio, -1.09))));
	}

	public static String ProbFailure(final double disp)
	{
		DecimalFormat fmt = new DecimalFormat(fmtThree);
		return fmt.format(0.335 * (1 - Math.exp(-0.048 * Math.pow(disp, 1.565))));
	}

	public static String[] BrayAndRathje(final double ky, final double h, final double vs, final double m, final double rock, final double r, final double mheaS, final double meanperS, final double sigdurS, final double normdispS, final double allowdisp, final boolean doScreening)
	{
		String ret[] = new String[13];

		double siteper, nrffact, meanper, dur, tstm, mheamhanrf, kmax, kykmax, normdisp, dispcm, dispin;
		double dur1, dur2, dur3, dur4, dur5;
		double arg;

		siteper = 4.0 * h / vs;
		nrffact = 0.62247 + 0.91958 * Math.exp(-rock / 0.44491);

		if(m <= 7.25)
		{
			meanper = (0.411 + 0.0837 * (m - 6.0) + 0.00208 * r) * Math.exp(meanperS * 0.437);
		}
		else
		{
			meanper = (0.411 + 1.25 * 0.0837 + 0.00208 * r) * Math.exp(meanperS * 0.437);
		}

		dur1 = Math.exp(-0.532 + 0.552 * Math.log((0.95 - 0.05) / (1.0 - 0.95)) - 0.0262 * Math.pow(Math.log((0.95 - 0.05) / (1 - 0.95)), 2));
		dur2 = Math.exp(5.204 + 0.851 * (m - 6.0));
		dur3 = Math.pow(10, (1.5 * m + 16.05));
		dur4 = Math.pow(dur2 / dur3, -1.0 / 3.0) / (4900000 * 3.2) + 0.063 * (r - 10.0);
		dur5 = Math.pow(dur2 / dur3, -1.0 / 3.0) / (4900000 * 3.2);

		if(r >= 10.0)
		{
			dur = Math.exp(Math.log(Math.exp(Math.log(dur4) + Math.log(dur1) + (0.493 * sigdurS))));
		}
		else
		{
			dur = Math.exp(Math.log(Math.exp(Math.log(dur5) + Math.log(dur1) + (0.493 * sigdurS))));
		}

		arg = (siteper / meanper);
		tstm = arg > 8.0 ? 8.0 : arg;

		arg = Math.exp(-0.6244 - 0.7831 * Math.log(tstm) + 0.298 * mheaS);
		mheamhanrf = arg > 1.0 ? 1.0 : arg;

		kmax = mheamhanrf * rock * nrffact;

		kykmax = ky / kmax;

		normdisp = Math.pow(10, 1.87 - 3.477 * kykmax + (normdispS * 0.35));

		dispcm = normdisp * dur * kmax;
		dispin = dispcm / 2.54;

		int incr = 0;
		DecimalFormat fmt = new DecimalFormat(fmtThree);
		DecimalFormat fmt1 = new DecimalFormat(fmtOne);
		ret[incr++] = fmt.format(siteper);
		ret[incr++] = fmt.format(nrffact);
		ret[incr++] = fmt.format(meanper);
		ret[incr++] = fmt.format(dur);
		ret[incr++] = fmt.format(tstm);
		ret[incr++] = fmt.format(mheamhanrf);
		ret[incr++] = fmt.format(kmax);
		ret[incr++] = fmt.format(kykmax);
		ret[incr++] = fmt.format(normdisp);
		ret[incr++] = fmt1.format(dispcm);
		ret[incr++] = fmt1.format(dispin);

		if(doScreening)
		{
			double medianfreq = nrffact / 3.477 * (1.87 - log10(allowdisp / (rock * nrffact * dur)));
			ret[incr++] = fmt.format(medianfreq); // medianfreq
			ret[incr++] = fmt.format(medianfreq * rock); // seiscoef
		}
		else
		{
			ret[incr++] = "";
			ret[incr++] = "";
		}

		return ret;
	}

	public static String Redigitize(DoubleList data, FileWriter ofile, final double di) throws IOException
	{
		Double val;
		double d, r, u, t1 = 0, t2, a1 = 0, a2, t0, a0;
		d = di;
		data.reset();
		if((val = data.each()) == null)
		{
			ofile.close();
			return "No data";
		}
		t2 = val.doubleValue();
		if((val = data.each()) == null)
		{
			ofile.close();
			return "Odd number of values";
		}
		a2 = val.doubleValue();
		boolean flag = false;
		for(int i = 0; !flag; i++)
		{
			t0 = (double)i * d;
			while(t0 > t2)
			{
				t1 = t2;
				a1 = a2;
				if((val = data.each()) == null)
				{
					flag = true;
					break;
				}
				t2 = val.doubleValue();
				if((val = data.each()) == null)
				{
					flag = true;
					return "Odd number of values";
				}
				a2 = val.doubleValue();
			}
			if(flag) break;
			r = t2 - t1;
			u = a2 - a1;
			if(r == 0)
			{
				a0 = a1;
			}
			else
			{
				a0 = a1 + (t0 - t1) * u / r;
			}
			ofile.write(Double.toString(a0));
			ofile.write('\n');
		}
		ofile.close();
		return "";
	}

	public static String TotalDuration(DoubleList data, final double di)
	{
		DecimalFormat fmt = new DecimalFormat(fmtOne);
		return fmt.format((double)(data.size()) * di);
	}

	public static String MeanPer(DoubleList data, final double di)
	{
		double[] arr = new double[data.size()];

		Double temp;
		data.reset();
		for(int i = 0; (temp = data.each()) != null; i++)
			arr[i] = temp.doubleValue();

		rdc(arr);

		taper(arr);

		// Pads the array so its length is a power of 2.
		int test = 0;

		for(int i = 1; test < arr.length; i++)
		{
			test = (int)Math.pow(2, i);
		}

		double[][] narr = new double[test][2];

		for(int i = 0; i < arr.length; i++)
		{
			narr[i][0] = arr[i];
			narr[i][1] = 0;
		}

		for(int i = narr.length; i < test; i++)
		{
			narr[i][0] = 0;
			narr[i][1] = 0;
		}

		// forward fft
		fft(narr);

		// scale to keep units correct
		for(int i = 0; i < narr.length; i++)
		{
			narr[i][0] *= di;
			narr[i][1] *= di;
		}

		// set frequency increment
		double df = 1.0 / ((double)(narr.length) * di);

		double top = 0, bot = 0, top2 = 0, bot2 = 0, tms = 0, to = 0;
		double f;

		for(int i = 0; i < arr.length; i++)
		{
			arr[i] = Math.sqrt(Math.pow(narr[i][0], 2) + Math.pow(narr[i][1], 2));
			f = i * df;

			if(f > 0.25 && f < 20.0)
			{
				top += (1.0 / f) * Math.pow(arr[i], 2);
				bot += Math.pow(arr[i], 2);
				top2 += Math.pow(1.0 / f, 2) * Math.pow(arr[i], 2);
				bot2 += Math.pow(arr[i], 2);
			}
		}

		DecimalFormat fmt = new DecimalFormat(fmtTwo);
		return fmt.format(top / bot);
	}

	public static void fft1(double[][] arr)
	{
		double temp[], carg, cw;
		temp = new double[2];
		int lx, m, j, l, istep;

		lx = arr.length;
		j = 1;

		for(int i = 1; i < lx; i++)
		{
			if(i < j)
			{
				temp[0] = arr[j - 1][0];
				temp[1] = arr[j - 1][1];

				arr[j - 1][0] = arr[i - 1][0];
				arr[j - 1][1] = arr[i - 1][1];

				arr[i - 1][0] = temp[0];
				arr[i - 1][1] = temp[1];
			}

			m = lx / 2;

			while(m < j)
			{
				j -= m;
				m /= 2;
			}

			j += m;
		}

		l = 1;

		do
		{
			istep = l + 1;

			for(m = 0; m < l; m++)
			{
				carg = (-Math.PI) * (double)(m - 1) / (double)l;
				cw = Math.exp(carg);
				for(int i = m; i < lx; lx += istep)
				{
					temp[0] = arr[i + 1][0];
					temp[1] = cw * arr[i + 1][1];

					arr[i + 1][0] = arr[i][0] - temp[0];
					arr[i + 1][1] = arr[i][1] - temp[1];

					arr[i][0] += temp[0];
					arr[i][1] += temp[1];
				}
			}

			l = istep;
		}	while(l < lx);
	}

	public static void fft(double[][] array)
	{
		double u_r,u_i, w_r,w_i, t_r,t_i;
		int ln, nv2, k, l, le, le1, j, ip, i, n;

		n = array.length;
		ln = (int)(Math.log((double)n) / Math.log(2) + 0.5);
		nv2 = n / 2;
		j = 1;

		for (i = 1; i < n; i++ )
		{
			if (i < j)
			{
				t_r = array[i - 1][0];
				t_i = array[i - 1][1];
				array[i - 1][0] = array[j - 1][0];
				array[i - 1][1] = array[j - 1][1];
				array[j - 1][0] = t_r;
				array[j - 1][1] = t_i;
			}

			k = nv2;

			while (k < j)
			{
				j = j - k;
				k = k / 2;
			}

			j = j + k;
		}

		for (l = 1; l <= ln; l++) /* loops thru stages */
		{
			le = (int)(Math.exp((double)l * Math.log(2)) + 0.5);
			le1 = le / 2;
			u_r = 1.0;
			u_i = 0.0;
			w_r =  Math.cos(Math.PI / (double)le1);
			w_i = -Math.sin(Math.PI / (double)le1);

			for (j = 1; j <= le1; j++) /* loops thru 1/2 twiddle values per stage */
			{
				for (i = j; i <= n; i += le) /* loops thru points per 1/2 twiddle */
				{
					ip = i + le1;
					t_r = array[ip - 1][0] * u_r - u_i * array[ip - 1][1];
					t_i = array[ip - 1][1] * u_r + u_i * array[ip - 1][0];

					array[ip - 1][0] = array[i - 1][0] - t_r;
					array[ip - 1][1] = array[i - 1][1] - t_i;

					array[i - 1][0] =  array[i - 1][0] + t_r;
					array[i - 1][1] =  array[i - 1][1] + t_i;
				}
				t_r = u_r * w_r - w_i * u_i;
				u_i = w_r * u_i + w_i * u_r;
				u_r = t_r;
			}
		}
	}


	// Removes a dc shift from the data by removing the mean.
	public static void rdc(double[] arr)
	{
		double sum = 0, mean;

		for(int i = 0; i < arr.length; i++)
			sum += arr[i];

		mean = sum / (double)(arr.length);

		for(int i = 0; i < arr.length; i++)
			arr[i] -= mean;
	}

	// Tapers the first and last 5% of the data.
	public static void taper(double[] arr)
	{
		double n, arg;
		double taper = 5;

		n = (double)(arr.length) * taper / 100.0;

		// beginning 5%
		for(int i = 0; i < n; i++)
		{
			arg = Math.PI * (double)(i) / (double)(n) + Math.PI;
			arr[i] *= (1.0 + Math.cos(arg)) / 2;
		}

		//ending 5%
		for(int i = 0; i < n; i++)
		{
			arg = Math.PI * (double)(i) / (double)(n) + Math.PI;
			arr[arr.length - i - 1] *= (1.0 + Math.cos(arg)) / 2;
		}
	}

	// coupled

	public static String Coupled(final DoubleList data, final double G, final double di, final double scale, final double uwgt, final double height, final double vs, final double damp, final double angle, final double caList[][])
	{
		double temp[];
		final double beta = 0.25;
		final double gamma = 0.5;
		double Mtot, M, L, omega;
		int qq;

		double disp[] = new double[caList.length];
		double mu[] = new double[caList.length];

		for(int ii = 0; ii < caList.length; ii++)
		{
			disp[ii] = caList[ii][0];
			mu[ii] = caList[ii][1];
		}

		double rho, time = 0;
		int i, k, j, kk;
		boolean slide;

		final double angleR = Math.toRadians(angle);
		final double angleC = Math.cos(angleR);
		final double angleS = Math.sin(angleR);

		setValueSize(di);

		data.reset();

		Double val;
		double cur = 0, prev;
		double curG = 0, prevG;

		// slide=0 no sliding, slide=1 sliding
		// variables that end in 1 are for previous time step
		// variables that end in 2 are for current time step

		double s1=0, sdot1=0, sdotdot1=0;
		double s2=0, sdot2=0, sdotdot2=0;
		double u1=0, udot1=0, udotdot1=0;
		double u2=0, udot2=0, udotdot2=0, baseacc=0;
		double basef=0, acc1=0, acc2=0, normalf1=0, normalf2=0;

		rho = uwgt / G;

		// calculate constants for Newmark algorithm

		Mtot = rho * height;
		slide = false;

		// qq indicates which mu is in effect
		qq = 0;

		omega = Math.PI * vs / (2.0 * height);
		L = 2.0 * rho * height / Math.PI;
		M = rho * height / 2.0;

		for(j = 0; (val = data.each()) != null; j++)
		{
			prev = cur;
			prevG = curG;
			cur = val.doubleValue() * scale;
			curG = cur * G;

			// setup state frm previous time step
			if(j == 0) // first time only
			{
				u1 = 0;
				udot1 = 0;
				udotdot1 = 0;
				s1 = 0;
				sdot1 = 0;
				sdotdot1 = 0;
				normalf1 = 0;
			}
			else
			{
				u1 = u2;
				udot1 = udot2;
				udotdot1 = udotdot2;
				s1 = s2;
				sdot1 = sdot2;
				sdotdot1 = sdotdot2;
				normalf1 = normalf2;
			}

			// setup acceleration loading

			// normal force corrected for vertical component of accel
			normalf2 = Mtot * G * angleC + Mtot * curG * angleS;

			if(j == 0)
			{
				acc1 = 0;
				acc2 = curG * angleC;
			}
			else
			{
				if(!slide)
				{
					acc1 = prevG * angleC;
					acc2 = curG * angleC;
				}
				else
				{
					acc1 = G * angleS - mu[qq] * normalf1 / Mtot;
					acc2 = G * angleS - mu[qq] * normalf2 / Mtot;
				}
			}

			// solve for u, udot, udotdot at next time step

			temp = CoupledSolvu(u1, udot1, udotdot1, u2, udot2, udotdot2, acc1, acc2, slide, j, M, Mtot, L, omega, beta, gamma, di, damp, G);

			u1 = temp[0];
			udot1 = temp[1];
			udotdot1 = temp[2];
			u2 = temp[3];
			udot2 = temp[4];
			udotdot2 = temp[5];
			acc1 = temp[6];
			acc2 = temp[7];

			// calculate base force based on udotdot calculation

			basef = -Mtot * curG * angleC - L * udotdot2 + Mtot * G * angleS;

			// check if sliding has started

			if(!slide)
			{
				if(basef > mu[qq] * normalf2)
					slide = true;
			}


			// based oncalculated response:

			if(slide)
			{
				// update sliding acceleration
				sdotdot2 = -curG * angleC - mu[qq] * normalf2 / Mtot - L * udotdot2 / Mtot + G * angleS;

				// if sliding is occuring, integrate sdotdot using trapezoid rule to get sdot and s
				sdot2 = sdot1 + 0.5 * di * (sdotdot2 + sdotdot1);
				s2 = s1 + 0.5 * di * (sdot2 + sdot1);

				// check if sliding has stopped
				if(sdot2 <= 0.0)
				{
					temp = CoupledSlideStop(s1, sdot1, sdotdot1, sdotdot2, u1, udot1, udotdot1, s2, sdot2, u2, udot2, udotdot2, slide, normalf2, Mtot, M, j, L, omega, mu[qq], beta, gamma, di, curG, prevG, angleS, angleC, damp, G);
					s1 = temp[0];
					sdot1 = temp[1];
					sdotdot1 = temp[2];
					u1 = temp[3];
					udot1 = temp[4];
					udotdot1 = temp[5];
					s2 = temp[6];
					sdot2 = temp[7];
					sdotdot2 = temp[8];
					u2 = temp[9];
					udot2 = temp[10];
					udotdot2 = temp[11];
					normalf2 = temp[12];

					slide = false;
					sdot2 = 0;
					sdotdot2 = 0;
				}
			}

			baseacc = curG * angleC;

			// output sliding quantities

			store(s2);

			if(!slide && Math.abs(s2) >= disp[qq] && qq < (mu.length - 1))
			{
				qq++;
			}

			time += di;
		}

		end(s2);

		return Double.toString(s2);
	}

	public static double[] CoupledSolvu(double u1, double udot1, double udotdot1, double u2, double udot2, double udotdot2, double acc1, double acc2, boolean slide, final int j, final double M, final double Mtot, final double L, final double omega, final double beta, final double gamma, final double di, final double damp, final double G)
	{
		double khat, a, b, dip, diu, diudot;
		double d1;

		if(slide)
			d1 = 1.0 - (L * L) / (M * Mtot);
		else
			d1 = 1;

		khat = (omega * omega) + 2.0 * damp * omega * gamma / (beta * di) + d1 / (beta * (di * di));
		a = d1 / (beta * di) + 2 * damp * omega * gamma / beta;
		b = d1 / (2 * beta) + di * 2 * damp * omega * (gamma / (2 * beta) - 1);

		if(j == 0)
		{
			dip = -L / M * (acc2 - acc1);
			diu = dip / khat;
			diudot = gamma / (beta * di) * diu;
			u2 = diu;
			udot2 = diudot;
			udotdot2 = (-(L / M) * acc2 - 2.0 * damp * omega * udot2 - (omega * omega) * u2) / d1;
		}
		else
		{
			dip = -L / M * (acc2 - acc1) + a * udot1 + b * udotdot1;
			diu = dip / khat;
			diudot = gamma / (beta * di) * diu - gamma / beta * udot1 + di * (1.0 - gamma / (2.0 * beta)) * udotdot1;
			u2 = u1 + diu;
			udot2 = udot1 + diudot;
			udotdot2 = (-(L / M) * acc2 - 2.0 * damp * omega * udot2 - (omega * omega) * u2) / d1;
		}

		double ret[] = new double[8];
		ret[0] = u1;
		ret[1] = udot1;
		ret[2] = udotdot1;
		ret[3] = u2;
		ret[4] = udot2;
		ret[5] = udotdot2;
		ret[6] = acc1;
		ret[7] = acc2;

		return ret;
	}

	public static double[] CoupledSlideStop(double s1, double sdot1, double sdotdot1, double sdotdot2, double u1, double udot1, double udotdot1, double s2, double sdot2, double u2, double udot2, double udotdot2, boolean slide, double normalf2, final double Mtot, final double M, final int j, final double L, final double omega, final double mu, final double beta, final double gamma, final double di, final double curG, final double prevG, final double angleS, final double angleC, final double damp, final double G)
	{
		double ddt, acc1, acc2;
		double acc1b, dd;
		double khat, dip, a, b;

		// end of slide time is taken where sdot=0 from previous analysis
		// assumin sliding throughout the time step

		dd = -sdot1 / (sdot2 - sdot1);
		ddt = dd * di;
		acc1 = G * angleS - mu * (G * angleC + curG * angleS);
		acc1b = prevG + dd * (curG - prevG);
		acc2 = G * angleS - mu * (G * angleC + acc1b * angleS);


		 // if dd=0, sliding has already stopped: skip this solution

		if(dd !=  0)
		{
			slide  =  true;

			double temp[] = CoupledSolvu(u1, udot1, udotdot1, u2, udot2, udotdot2, acc1, acc2, slide, j, M, Mtot, L, omega, beta, gamma, di, damp, G);
			u1 = temp[0];
			udot1 = temp[1];
			udotdot1 = temp[2];
			u2 = temp[3];
			udot2 = temp[4];
			udotdot2 = temp[5];
			acc1 = temp[6];
			acc2 = temp[7];

			u1 = u2;
			udot1 = udot2;
			udotdot1 = udotdot2;
			normalf2 = Mtot * G * angleC + Mtot * acc1b * angleS;
			sdotdot2 =  - acc1b * angleC - mu * normalf2 / Mtot - L * udotdot2 / Mtot + G * angleS;
			sdot2 = sdot1 + 0.5 * ddt * (sdotdot2 + sdotdot1);
			s2 = s1 + 0.5 * ddt * (sdot1 + sdot2);

			 // solve for non sliding response during remaining part of di

			ddt = (1.0 - dd) * di;
			// slide = false; // this does nothing, afaik
			acc1 = acc2;
			acc2 = curG * angleC;

			khat = 1.0 + 2.0 * damp * omega * gamma * ddt + (omega * omega) * beta * (ddt * ddt);
			a = (1.0 - (L * L) / (Mtot * M)) + 2.0 * damp * omega * ddt * (gamma - 1.0) + (omega * omega) * (ddt * ddt) * (beta - 0.5);
			b = (omega * omega) * ddt;
			dip =  - L / M * (acc2 - acc1) + a * (udotdot1) - b * (udot1);
			udotdot2 = dip / khat;

			udot2 = udot1 + (1.0 - gamma) * ddt * (udotdot1) + gamma * ddt * (udotdot2);
			u2 = u1 + udot1 * ddt + (0.5 - beta) * (ddt * ddt) * (udotdot1) + beta * (ddt * ddt) * (udotdot2);
		}

		double ret[] = new double[13];
		ret[0] = s1;
		ret[1] = sdot1;
		ret[2] = sdotdot1;
		ret[3] = u1;
		ret[4] = udot1;
		ret[5] = udotdot1;
		ret[6] = s2;
		ret[7] = sdot2;
		ret[8] = sdotdot2;
		ret[9] = u2;
		ret[10] = udot2;
		ret[11] = udotdot2;
		ret[12] = normalf2;

		return ret;
	}
}
