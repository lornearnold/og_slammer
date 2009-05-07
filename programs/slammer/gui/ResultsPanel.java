/* This file is in the public domain. */

package slammer.gui;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.border.*;
import org.jfree.data.xy.*;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.chart.*;
import org.jfree.chart.axis.*;
import java.text.DecimalFormat;
import java.util.*;
import java.io.*;
import slammer.*;
import slammer.analysis.*;

class ResultsPanel extends JPanel implements ActionListener
{
	// array indexes
	public final static int RB = 0; // rigid block
	public final static int DC = 1; // decoupled
	public final static int CP = 2; // coupled
	public final static int ANALYSIS_TYPES = 3;
	public final static int AVG = 2; // average
	public final static int NOR = 0; // normal
	public final static int INV = 1; // inverse

	// table column indicies
	public final static int WIDTH = 3;
	public final static int OFFSET = 2;
	public final static int RBC = OFFSET;
	public final static int DCC = RBC + WIDTH;
	public final static int CPC = DCC + WIDTH;
	public final static int LEN = CPC + WIDTH;

	String polarityName[] = { "Normal", "Inverse", "Average" };

	SlammerTabbedPane parent;

	JTextField decimalsTF = new JTextField("1", 2);

	JButton Analyze = new JButton("Perform analyses");
	JButton ClearOutput = new JButton("Clear output");

	JEditorPane outputPane = new JEditorPane("text/html", "");
	JScrollPane outputPanes = new JScrollPane(outputPane);

	JButton saveResultsOutput = new JButton("Save results table");
	JFileChooser fc = new JFileChooser();

	JButton plotHistogram = new JButton("Plot histogram of displacements");
	JButton plotDisplacement = new JButton("Plot displacements versus time");
	JCheckBox plotDisplacementLegend = new JCheckBox("Display legend", false);
	JTextField outputBins = new JTextField("10", 2);

	JRadioButton polarityNorDisp = new JRadioButton(polarityName[NOR], true);
	JRadioButton polarityInvDisp = new JRadioButton(polarityName[INV]);
	ButtonGroup polarityGroupDisp = new ButtonGroup();

	JRadioButton polarityAvgHist = new JRadioButton(polarityName[AVG], true);
	JRadioButton polarityNorHist = new JRadioButton(polarityName[NOR]);
	JRadioButton polarityInvHist = new JRadioButton(polarityName[INV]);
	ButtonGroup polarityGroupHist = new ButtonGroup();

	JCheckBox analysisDisp[] = new JCheckBox[ANALYSIS_TYPES];
	JRadioButton analysisHist[] = new JRadioButton[ANALYSIS_TYPES];
	String analysisTitle[] = new String[ANALYSIS_TYPES];
	ButtonGroup analysisHistGroup = new ButtonGroup();

	ArrayList results;
	XYSeries xys[][][];

	JRadioButton outputDelTab = new JRadioButton("Text, tab delimited", true);
	JRadioButton outputDelSpace = new JRadioButton("Text, space delimited");
	JRadioButton outputDelComma = new JRadioButton("Text, comma delimited");
	JRadioButton outputDelHTML = new JRadioButton("HTML");
	ButtonGroup outputDelGroup = new ButtonGroup();

	boolean paramUnit = false;
	String unitDisplacement = "";
	DecimalFormat unitFmt;
	ArrayList dataVect[][];
	XYSeriesCollection xycol;
	String parameters;

	public ResultsPanel(SlammerTabbedPane parent) throws Exception
	{
		this.parent = parent;

		outputPane.setEditable(false);

		Analyze.setActionCommand("analyze");
		Analyze.addActionListener(this);

		ClearOutput.setActionCommand("clearOutput");
		ClearOutput.addActionListener(this);

		saveResultsOutput.setActionCommand("saveResultsOutput");
		saveResultsOutput.addActionListener(this);

		plotHistogram.setActionCommand("plotHistogram");
		plotHistogram.addActionListener(this);

		plotDisplacement.setActionCommand("plotDisplacement");
		plotDisplacement.addActionListener(this);

		plotDisplacementLegend.setActionCommand("plotDisplacementLegend");
		plotDisplacementLegend.addActionListener(this);

		outputDelGroup.add(outputDelTab);
		outputDelGroup.add(outputDelSpace);
		outputDelGroup.add(outputDelComma);
		outputDelGroup.add(outputDelHTML);

		polarityGroupDisp.add(polarityNorDisp);
		polarityGroupDisp.add(polarityInvDisp);

		polarityGroupHist.add(polarityAvgHist);
		polarityGroupHist.add(polarityNorHist);
		polarityGroupHist.add(polarityInvHist);

		analysisDisp[RB] = new JCheckBox(ParametersPanel.stringRB);
		analysisDisp[DC] = new JCheckBox(ParametersPanel.stringDC);
		analysisDisp[CP] = new JCheckBox(ParametersPanel.stringCP);
		analysisHist[RB] = new JRadioButton(ParametersPanel.stringRB);
		analysisHist[DC] = new JRadioButton(ParametersPanel.stringDC);
		analysisHist[CP] = new JRadioButton(ParametersPanel.stringCP);

		analysisHistGroup.add(analysisHist[RB]);
		analysisHistGroup.add(analysisHist[DC]);
		analysisHistGroup.add(analysisHist[CP]);

		analysisTitle[RB] = "Rigid-Block";
		analysisTitle[DC] = "Decoupled";
		analysisTitle[CP] = "Coupled";

		setLayout(new BorderLayout());

		add(BorderLayout.NORTH, createHeader());
		add(BorderLayout.CENTER, createTable());
		add(BorderLayout.SOUTH, createGraphs());

		clearOutput();
	}

	public void actionPerformed(java.awt.event.ActionEvent e)
	{
		try
		{
			String command = e.getActionCommand();
			if(command.equals("analyze"))
			{
				final SwingWorker worker = new SwingWorker()
				{
					ProgressFrame pm = new ProgressFrame(0);

					public Object construct()
					{
						try
						{
							clearOutput();

							paramUnit = parent.Parameters.unitMetric.isSelected();
							final double g = paramUnit ? Analysis.Gcmss : Analysis.Ginss;
							unitDisplacement = paramUnit ? "(cm)" : "(in.)";
							StringBuilder outputText = new StringBuilder();
							outputText.append("<html><body><table border=\"1\"><tr>" +
								"<th>Earthquake</th><th>Record</th>" +
								"<th colspan=\"3\">" + ParametersPanel.stringRB + " " + unitDisplacement + "</th>" +
								"<th colspan=\"3\">" + ParametersPanel.stringDC + " " + unitDisplacement + "</th>" +
								"<th colspan=\"3\">" + ParametersPanel.stringCP + " " + unitDisplacement + "</th>" +
								"</tr>");

							boolean paramDualslope = parent.Parameters.dualSlope.isSelected();
							Double d;

							double paramScale;
							if(parent.Parameters.scalePGA.isSelected())
							{
								d = (Double)Utils.checkNum(parent.Parameters.scalePGAval.getText(), "scale PGA field", null, false, null, new Double(0), false, null, false);
								if(d == null)
								{
									parent.selectParameters();
									return null;
								}
								paramScale = d.doubleValue();
							}
							else if(parent.Parameters.scaleOn.isSelected())
							{
								d = (Double)Utils.checkNum(parent.Parameters.scaleData.getText(), "scale data field", null, false, null, null, false, null, false);
								if(d == null)
								{
									parent.selectParameters();
									return null;
								}
								paramScale = d.doubleValue();
							}
							else
								paramScale = 0;

							changeDecimal();

							boolean paramRigid = parent.Parameters.typeRigid.isSelected();
							boolean paramDecoupled = parent.Parameters.typeDecoupled.isSelected();
							boolean paramCoupled = parent.Parameters.typeCoupled.isSelected();

							graphDisp(paramRigid, paramDecoupled, paramCoupled);

							if(!paramRigid && !paramDecoupled && !paramCoupled)
							{
								parent.selectParameters();
								GUIUtils.popupError("Error: no analyses methods selected.");
								return null;
							}

							Object[][] res = Utils.getDB().runQuery("select eq, record, digi_int, path, pga from data where select2=1 and analyze=1");

							if(res == null || res.length <= 1)
							{
								parent.selectSelectRecords();
								GUIUtils.popupError("No records selected for analysis.");
								return null;
							}

							xys = new XYSeries[res.length][3][2];
							dataVect = new ArrayList[3][3];

							String eq, record;
							DoubleList dat;
							double di;
							int num = 0;
							double avg;
							double total[][] = new double[3][3];
							double scale = 1, iscale, scaleRB;
							double inv, norm;
							double[][] ca;
							double[] ain = null;
							double thrust = 0, uwgt = 0, height = 0, vs = 0, damp = 0, vr = 0;
							boolean dv3 = false;

							scaleRB = paramUnit ? 1 : Analysis.CMtoIN;

							if(parent.Parameters.CAdisp.isSelected())
							{
								String value;
								java.util.Vector caVect;
								TableCellEditor editor = null;

								editor = parent.Parameters.dispTable.getCellEditor();
								caVect = parent.Parameters.dispTableModel.getDataVector();

								if(editor != null)
									editor.stopCellEditing();

								ca = new double[caVect.size()][2];

								for(int i = 0; i < caVect.size(); i++)
								{
									for(int j = 0; j < 2; j++)
									{
										value = (String)(((java.util.Vector)(caVect.get(i))).get(j));
										if(value == null || value == "")
										{
											parent.selectParameters();
											GUIUtils.popupError("Error: empty field in table.\nPlease complete the displacement table so that all data pairs have values, or delete all empty rows.");
											return null;
										}
										d = (Double)Utils.checkNum(value, "displacement table", null, false, null, null, false, null, false);
										if(d == null)
										{
											parent.selectParameters();
											return null;
										}
										ca[i][j] = d.doubleValue();
									}
								}

								if(caVect.size() == 0)
								{
									parent.selectParameters();
									GUIUtils.popupError("Error: no displacements listed in displacement table.");
									return null;
								}
							}
							else
							{
								d = (Double)Utils.checkNum(parent.Parameters.CAconstTF.getText(), "constant critical acceleration field", null, false, null, new Double(0), true, null, false);
								if(d == null)
								{
									parent.selectParameters();
									return null;
								}
								ca = new double[1][2];
								ca[0][0] = 0;
								ca[0][1] = d.doubleValue();
							}

							if(paramRigid && paramDualslope)
							{
								Double thrustD = (Double)Utils.checkNum(parent.Parameters.thrustAngle.getText(), "thrust angle field", new Double(90), true, null, new Double(0), true, null, false);
								if(thrustD == null)
								{
									parent.selectParameters();
									return null;
								}
								else
									thrust = thrustD.doubleValue();
							}

							if(paramDecoupled || paramCoupled)
							{
								Double tempd;

								uwgt = 100.0; // unit weight is disabled, so hardcode to 100

								tempd = (Double)Utils.checkNum(parent.Parameters.paramHeight.getText(), ParametersPanel.stringHeight + " field", null, false, null, null, false, null, false);
								if(tempd == null)
								{
									parent.selectParameters();
									return null;
								}
								else
									height = tempd.doubleValue();

								tempd = (Double)Utils.checkNum(parent.Parameters.paramVs.getText(), ParametersPanel.stringVs + " field", null, false, null, null, false, null, false);
								if(tempd == null)
								{
									parent.selectParameters();
									return null;
								}
								else
									vs = tempd.doubleValue();

								tempd = (Double)Utils.checkNum(parent.Parameters.paramDamp.getText(), ParametersPanel.stringDamp + " field", null, false, null, null, false, null, false);
								if(tempd == null)
								{
									parent.selectParameters();
									return null;
								}
								else
									damp = tempd.doubleValue() / 100.0;

								dv3 = parent.Parameters.paramSoilModel.getSelectedIndex() == 1;

								// metric
								if(paramUnit)
								{
									uwgt /= Analysis.M3toCM3;
									height *= Analysis.MtoCM;
									vs *= Analysis.MtoCM;
									vr *= Analysis.MtoCM;
								}
								// english
								else
								{
									uwgt /= Analysis.FT3toIN3;
									height *= Analysis.FTtoIN;
									vs *= Analysis.FTtoIN;
									vr *= Analysis.FTtoIN;
								}
							}

							File testFile;
							String path;

							if(paramRigid)
							{
								dataVect[RB][NOR] = new ArrayList<Double>(res.length - 1);
								dataVect[RB][INV] = new ArrayList<Double>(res.length - 1);
								dataVect[RB][AVG] = new ArrayList<Double>(res.length - 1);
							}

							if(paramDecoupled)
							{
								dataVect[DC][NOR] = new ArrayList<Double>(res.length - 1);
								dataVect[DC][INV] = new ArrayList<Double>(res.length - 1);
								dataVect[DC][AVG] = new ArrayList<Double>(res.length - 1);
							}

							if(paramCoupled)
							{
								dataVect[CP][NOR] = new ArrayList<Double>(res.length - 1);
								dataVect[CP][INV] = new ArrayList<Double>(res.length - 1);
								dataVect[CP][AVG] = new ArrayList<Double>(res.length - 1);
							}

							iscale = -1.0 * scale;

							pm.setMaximum(res.length);

							int j, k;
							Object[] row;
							results = new ArrayList();

							results.add(new Object[] {
								"Earthquake", "Record",
								"<---", ParametersPanel.stringRB + " " + unitDisplacement, "---->",
								"<---", ParametersPanel.stringDC + " " + unitDisplacement, "---->",
								"<---", ParametersPanel.stringCP + " " + unitDisplacement, "---->"
							});

							row = new Object[] { null, "Polarity:",
								polarityName[NOR], polarityName[INV], polarityName[AVG],
								polarityName[NOR], polarityName[INV], polarityName[AVG],
								polarityName[NOR], polarityName[INV], polarityName[AVG] };

							outputText.append(makeRow(row, true));
							results.add(row);

							for(int i = 1; i < res.length && !pm.isCanceled(); i++)
							{
								row = new Object[LEN];
								eq = res[i][0].toString();
								record = res[i][1].toString();

								pm.update(i, eq + " - " + record);

								row[0] = eq;
								row[1] = record;

								path = res[i][3].toString();
								testFile = new File(path);
								if(!testFile.exists() || !testFile.canRead())
								{
									row[2] = "File does not exist or is not readable";
									row[3] = path;
									outputText.append(makeRow(row));
									results.add(row);
									continue;
								}

								dat = new DoubleList(path, 0, parent.Parameters.scaleOn.isSelected() ? paramScale : 1.0);
								if(dat.bad())
								{
									row[2] = "Invalid data at point " + dat.badEntry();
									row[3] = path;
									outputText.append(makeRow(row));
									results.add(row);
									continue;
								}

								num++;

								di = Double.parseDouble(res[i][2].toString());

								if(parent.Parameters.scalePGA.isSelected())
								{
									scale = paramScale / Double.parseDouble(res[i][4].toString());
									iscale = -scale;
								}

								if(paramCoupled || paramDecoupled)
									ain = dat.getAsArray();

								// do the actual analysis

								if(paramRigid)
								{
									norm = RigidBlock.SlammerRigorous("norm", dat, di, ca, scale, paramDualslope, thrust, scaleRB);
									Analysis.graphData.setKey(row[0] + " - " + row[1] + " - " + ParametersPanel.stringRB + ", " + polarityName[NOR]);
									xys[i - 1][RB][NOR] = Analysis.graphData;

									inv = RigidBlock.SlammerRigorous("inv", dat, di, ca, iscale, paramDualslope, thrust, scaleRB);
									Analysis.graphData.setKey(row[0] + " - " + row[1] + " - " + ParametersPanel.stringRB + ", " + polarityName[INV]);
									xys[i - 1][RB][INV] = Analysis.graphData;

									avg = avg(inv, norm);

									total[RB][NOR] += norm;
									total[RB][INV] += inv;
									total[RB][AVG] += avg;

									for(j = 0; j < dataVect[RB][NOR].size() && ((Double)dataVect[RB][NOR].get(j)).doubleValue() < norm; j++)
										;
									dataVect[RB][NOR].add(j, new Double(norm));

									for(j = 0; j < dataVect[RB][INV].size() && ((Double)dataVect[RB][INV].get(j)).doubleValue() < inv; j++)
										;
									dataVect[RB][INV].add(j, new Double(inv));

									for(j = 0; j < dataVect[RB][AVG].size() && ((Double)dataVect[RB][AVG].get(j)).doubleValue() < avg; j++)
										;
									dataVect[RB][AVG].add(j, new Double(avg));

									row[RBC + AVG] = unitFmt.format(avg);
									row[RBC + NOR] = unitFmt.format(norm);
									row[RBC + INV] = unitFmt.format(inv);
								}

								if(paramDecoupled)
								{
									// [i]scale is divided by Gcmss because the algorithm expects input data in Gs, but our input files are in cmss. this has nothing to do with, and is not affected by, the unit base being used (english or metric).
									norm = Decoupled.Decoupled(ain, uwgt, height, vs, damp, di, scale / Analysis.Gcmss, g, vr, ca, dv3);
									Analysis.graphData.setKey(row[0] + " - " + row[1] + " - " + ParametersPanel.stringDC + ", " + polarityName[NOR]);
									xys[i - 1][DC][NOR] = Analysis.graphData;

									inv = Decoupled.Decoupled(ain, uwgt, height, vs, damp, di, iscale / Analysis.Gcmss, g, vr, ca, dv3);
									Analysis.graphData.setKey(row[0] + " - " + row[1] + " - " + ParametersPanel.stringDC + ", " + polarityName[INV]);
									xys[i - 1][DC][INV] = Analysis.graphData;

									avg = avg(inv, norm);

									total[DC][NOR] += norm;
									total[DC][INV] += inv;
									total[DC][AVG] += avg;

									for(j = 0; j < dataVect[DC][AVG].size() && ((Double)dataVect[DC][AVG].get(j)).doubleValue() < avg; j++)
										;

									dataVect[DC][AVG].add(j, new Double(avg));
									dataVect[DC][NOR].add(j, new Double(norm));
									dataVect[DC][INV].add(j, new Double(inv));

									row[DCC + AVG] = unitFmt.format(avg);
									row[DCC + NOR] = unitFmt.format(norm);
									row[DCC + INV] = unitFmt.format(inv);
								}

								if(paramCoupled)
								{
									// [i]scale is divided by Gcmss because the algorithm expects input data in Gs, but our input files are in cmss. this has nothing to do with, and is not affected by, the unit base being used (english or metric).
									norm = Coupled.Coupled(ain, uwgt, height, vs, damp, di, scale / Analysis.Gcmss, g, vr, ca, dv3);
									Analysis.graphData.setKey(row[0] + " - " + row[1] + " - " + ParametersPanel.stringCP + ", " + polarityName[NOR]);
									xys[i - 1][CP][NOR] = Analysis.graphData;

									inv = Coupled.Coupled(ain, uwgt, height, vs, damp, di, iscale / Analysis.Gcmss, g, vr, ca, dv3);
									Analysis.graphData.setKey(row[0] + " - " + row[1] + " - " + ParametersPanel.stringCP + ", " + polarityName[INV]);
									xys[i - 1][CP][INV] = Analysis.graphData;

									avg = avg(inv, norm);

									total[CP][NOR] += norm;
									total[CP][INV] += inv;
									total[CP][AVG] += avg;

									for(j = 0; j < dataVect[CP][AVG].size() && ((Double)dataVect[CP][AVG].get(j)).doubleValue() < avg; j++)
										;

									dataVect[CP][AVG].add(j, new Double(avg));
									dataVect[CP][NOR].add(j, new Double(norm));
									dataVect[CP][INV].add(j, new Double(inv));

									row[CPC + AVG] = unitFmt.format(avg);
									row[CPC + NOR] = unitFmt.format(norm);
									row[CPC + INV] = unitFmt.format(inv);
								}

								outputText.append(makeRow(row));
								results.add(row);
							}
							pm.update("Calculating stastistics...");

							double mean, value, valtemp;
							Object[] rmean = new Object[LEN];
							Object[] rmedian = new Object[LEN];
							Object[] rsd = new Object[LEN];

							rmean[1] = "Mean value";
							rmedian[1] = "Median value";
							rsd[1] = "Standard deviation";

							for(j = 0; j < total.length; j++)
							{
								for(k = 0; k < total[j].length; k++)
								{
									if(dataVect[j][k] == null || dataVect[j][k].size() == 0)
										continue;

									mean = Double.parseDouble(unitFmt.format(total[j][k] / num));
									rmean[j * WIDTH + OFFSET + k] = unitFmt.format(mean);

									if(num % 2 == 0)
									{
										double fst = (Double)dataVect[j][k].get(num / 2);
										double snd = (Double)dataVect[j][k].get(num / 2 - 1);
										rmedian[j * WIDTH + OFFSET + k] = unitFmt.format(avg(fst, snd));
									}
									else
										rmedian[j * WIDTH + OFFSET + k] = unitFmt.format(dataVect[j][k].get(num / 2));

									value = 0;

									for(int i = 0; i < num; i++)
									{
										valtemp = mean - ((Double)dataVect[j][k].get(i)).doubleValue();
										value += (valtemp * valtemp);
									}

									value /= num;
									value = Math.sqrt(value);
									rsd[j * WIDTH + OFFSET + k] = unitFmt.format(value);
								}
							}

							outputText.append(makeRow(new Object[LEN]));
							outputText.append(makeRow(rmean));
							outputText.append(makeRow(rmedian));
							outputText.append(makeRow(rsd));
							outputText.append("</table></body></html>");

							results.add(new Object[LEN]);
							results.add(rmean);
							results.add(rmedian);
							results.add(rsd);

							outputPane.setText(outputText.toString());
						}
						catch(Throwable ex)
						{
							Utils.catchException(ex);
						}

						return null;
					}

					public void finished() {
						pm.dispose();
					}
				};
				worker.start();
			}
			else if(command.equals("clearOutput"))
			{
				clearOutput();
			}
			else if(command.equals("saveResultsOutput"))
			{
				if(fc.showSaveDialog(this) == JFileChooser.APPROVE_OPTION)
				{
					FileWriter fw = new FileWriter(fc.getSelectedFile());

					if(outputDelHTML.isSelected())
					{
						fw.write(outputPane.getText());
						fw.close();
					}
					else
					{
						String delim;
						if(outputDelSpace.isSelected())
							delim = " ";
						else if(outputDelComma.isSelected())
							delim = ",";
						else
							delim = "\t";

						Object o;
						Object[] r;

						for(int i = 0; i < results.size(); i++)
						{
							r = (Object[])results.get(i);

							for(int j = 0; j < r.length; j++)
							{
								if(j != 0)
									fw.write(delim);

								o = r[j];
								if(o == null)
									o = "";

								fw.write(o.toString());
							}

							fw.write("\n");
						}

						fw.close();
					}
				}
			}
			else if(command.equals("plotHistogram"))
			{
				if(dataVect == null)
					return;

				String name = "", title, pname;
				HistogramDataset dataset = new HistogramDataset();

				int polarity, analysis = -1;

				if(polarityAvgHist.isSelected()) polarity = AVG;
				else if(polarityNorHist.isSelected()) polarity = NOR;
				else if(polarityInvHist.isSelected()) polarity = INV;
				else polarity = -1;

				for(int i = 0; i < analysisHist.length; i++)
					if(analysisHist[i].isSelected())
						analysis = i;

				if(analysis == -1)
					return;

				pname = polarityName[polarity];

				Double Bins = (Double)Utils.checkNum(outputBins.getText(), "output bins field", null, false, null, new Double(0), false, null, false);

				if(Bins == null || dataVect[analysis][polarity] == null)
					return;

				name = analysisTitle[analysis];

				double series[] = new double[dataVect[analysis][polarity].size()];

				for(int j = 0; j < dataVect[analysis][polarity].size(); j++)
					series[j] = (((Double)dataVect[analysis][polarity].get(j)).doubleValue());

				dataset.addSeries(name, series, (int)Bins.doubleValue());

				title = "Histogram of " + name + " Displacements (" + pname + " Polarity)";

				JFreeChart hist = ChartFactory.createHistogram(title, "Displacement " + unitDisplacement, "Number of Records", dataset, org.jfree.chart.plot.PlotOrientation.VERTICAL, false, true, false);
				ChartFrame frame = new ChartFrame(title, hist);

				frame.pack();
				frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
			else if(command.equals("plotDisplacement"))
			{
				XYSeriesCollection xysc = new XYSeriesCollection();

				int polarity = polarityNorDisp.isSelected() ? NOR : INV;
				String pname = polarityName[polarity];

				String name = "";
				boolean first = true;

				for(int i = 0; i < analysisDisp.length; i++)
				{
					if(analysisDisp[i].isSelected() && dataVect[i][polarity] != null)
					{
						if(first)
							first = false;
						else
							name += ", ";

						name += analysisTitle[i];

						for(int j = 0; j < dataVect[i][polarity].size(); j++)
							xysc.addSeries(xys[j][i][polarity]);
					}
				}

				if(first)
					return;

				name += " Displacement versus Time";

				JFreeChart chart = ChartFactory.createXYLineChart(name, "Time (s)", "Displacement--" + pname + " polarity " + unitDisplacement, xysc, org.jfree.chart.plot.PlotOrientation.VERTICAL, plotDisplacementLegend.isSelected(), true, false);

				//margins are stupid
				chart.getXYPlot().getDomainAxis().setLowerMargin(0);
				chart.getXYPlot().getDomainAxis().setUpperMargin(0);
				chart.getXYPlot().getDomainAxis().setLowerBound(0);

				ChartFrame frame = new ChartFrame(name, chart);
				frame.pack();
				frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
		}
		catch (Exception ex)
		{
			Utils.catchException(ex);
		}
	}

	private JPanel createHeader()
	{
		JPanel containerL = new JPanel();
		containerL.add(new JLabel("Display results to"));
		containerL.add(decimalsTF);
		containerL.add(new JLabel("decimals. "));

		JPanel container = new JPanel(new BorderLayout());
		container.add(containerL, BorderLayout.WEST);
		container.add(Analyze, BorderLayout.CENTER);
		container.add(ClearOutput, BorderLayout.EAST);

		return container;
	}

	private JPanel createTable()
	{
		JPanel outputTablePanel = new JPanel(new BorderLayout());

		outputTablePanel.add(BorderLayout.CENTER, outputPanes);

		return outputTablePanel;
	}

	private JTabbedPane createGraphs()
	{
		JTabbedPane jtp = new JTabbedPane();

		jtp.add("Graph displacements", createGraphDisplacementsPanel());
		jtp.add("Graph histogram", createGraphHistogramPanel());
		jtp.add("Save results table", createSaveResultsPanel());

		return jtp;
	}

	private JPanel createGraphDisplacementsPanel()
	{
		JPanel panel = new JPanel();

		GridBagLayout gridbag = new GridBagLayout();
		panel.setLayout(gridbag);
		GridBagConstraints c = new GridBagConstraints();
		JLabel label;

		int x = 0;
		int y = 0;

		c.anchor = GridBagConstraints.NORTHWEST;

		c.gridx = x++;
		c.gridy = y++;
		label = new JLabel("Analyses:");
		gridbag.setConstraints(label, c);
		panel.add(label);

		c.gridy = y++;
		gridbag.setConstraints(analysisDisp[RB], c);
		panel.add(analysisDisp[RB]);

		c.gridy = y++;
		gridbag.setConstraints(analysisDisp[DC], c);
		panel.add(analysisDisp[DC]);

		c.gridy = y++;
		gridbag.setConstraints(analysisDisp[CP], c);
		panel.add(analysisDisp[CP]);

		y = 0;

		c.gridy = y++;
		c.gridx = x++;
		c.insets = new Insets(0, 5, 0, 0);
		label = new JLabel("Polarity:");
		gridbag.setConstraints(label, c);
		panel.add(label);

		c.gridy = y++;
		gridbag.setConstraints(polarityNorDisp, c);
		panel.add(polarityNorDisp);

		c.gridy = y++;
		gridbag.setConstraints(polarityInvDisp, c);
		panel.add(polarityInvDisp);

		y = 0;

		c.gridy = y++;
		c.gridx = x++;
		c.gridheight = 2;
		c.anchor = GridBagConstraints.WEST;
		gridbag.setConstraints(plotDisplacement, c);
		panel.add(plotDisplacement);

		y++;
		c.gridy = y;
		gridbag.setConstraints(plotDisplacementLegend, c);
		panel.add(plotDisplacementLegend);

		c.gridx = x;
		label = new JLabel();
		c.weightx = 1;
		gridbag.setConstraints(label, c);
		panel.add(label);

		return panel;
	}

	private JPanel createGraphHistogramPanel()
	{
		JPanel panel = new JPanel();

		GridBagLayout gridbag = new GridBagLayout();
		panel.setLayout(gridbag);
		GridBagConstraints c = new GridBagConstraints();
		JLabel label;

		int x = 0;
		int y = 0;

		c.anchor = GridBagConstraints.NORTHWEST;

		c.gridx = x++;
		c.gridy = y++;
		label = new JLabel("Analysis:");
		gridbag.setConstraints(label, c);
		panel.add(label);

		c.gridy = y++;
		gridbag.setConstraints(analysisHist[RB], c);
		panel.add(analysisHist[RB]);

		c.gridy = y++;
		gridbag.setConstraints(analysisHist[DC], c);
		panel.add(analysisHist[DC]);

		c.gridy = y++;
		gridbag.setConstraints(analysisHist[CP], c);
		panel.add(analysisHist[CP]);

		y = 0;

		c.gridy = y++;
		c.gridx = x++;
		c.insets = new Insets(0, 5, 0, 0);
		label = new JLabel("Polarity:");
		gridbag.setConstraints(label, c);
		panel.add(label);

		c.gridy = y++;
		gridbag.setConstraints(polarityAvgHist, c);
		panel.add(polarityAvgHist);

		c.gridy = y++;
		gridbag.setConstraints(polarityNorHist, c);
		panel.add(polarityNorHist);

		c.gridy = y++;
		gridbag.setConstraints(polarityInvHist, c);
		panel.add(polarityInvHist);

		y = 0;

		c.gridx = x++;
		c.gridy = y++;
		c.gridheight = 2;
		c.anchor = GridBagConstraints.WEST;
		gridbag.setConstraints(plotHistogram, c);
		panel.add(plotHistogram);

		JPanel bins = new JPanel();
		bins.add(new JLabel("Plot with "));
		bins.add(outputBins);
		bins.add(new JLabel(" bins."));
		y++;
		c.gridy = y;
		gridbag.setConstraints(bins, c);
		panel.add(bins);

		c.gridx = x;
		label = new JLabel();
		c.weightx = 1;
		gridbag.setConstraints(label, c);
		panel.add(label);

		return panel;
	}

	private JPanel createSaveResultsPanel()
	{
		JPanel panel = new JPanel();

		GridBagLayout gridbag = new GridBagLayout();
		panel.setLayout(gridbag);
		GridBagConstraints c = new GridBagConstraints();
		JLabel label;

		int x = 0;
		int y = 0;

		c.anchor = GridBagConstraints.WEST;

		c.gridx = x++;
		c.gridy = y++;
		label = new JLabel("Output format:");
		gridbag.setConstraints(label, c);
		panel.add(label);

		c.gridy = y++;
		gridbag.setConstraints(outputDelTab, c);
		panel.add(outputDelTab);

		c.gridy = y++;
		gridbag.setConstraints(outputDelSpace, c);
		panel.add(outputDelSpace);

		c.gridy = y++;
		gridbag.setConstraints(outputDelComma, c);
		panel.add(outputDelComma);

		c.gridy = y++;
		gridbag.setConstraints(outputDelHTML, c);
		panel.add(outputDelHTML);

		c.gridx = x++;
		c.gridheight = y;
		c.gridy = 0;
		gridbag.setConstraints(saveResultsOutput, c);
		panel.add(saveResultsOutput);

		c.gridx = x;
		label = new JLabel();
		c.weightx = 1;
		gridbag.setConstraints(label, c);
		panel.add(label);

		return panel;
	}

	private void clearOutput()
	{
		dataVect = null;
		xys = null;
		graphDisp(false, false, false);
		outputPane.setText(null);
	}

	private void graphDisp(boolean rigid, boolean decoupled, boolean coupled)
	{
		boolean any = rigid || decoupled || coupled;

		analysisDisp[RB].setEnabled(rigid);
		analysisHist[RB].setEnabled(rigid);
		analysisDisp[DC].setEnabled(decoupled);
		analysisHist[DC].setEnabled(decoupled);
		analysisDisp[CP].setEnabled(coupled);
		analysisHist[CP].setEnabled(coupled);

		plotHistogram.setEnabled(any);
		plotDisplacement.setEnabled(any);
	}

	private void changeDecimal()
	{
		int i;

		Double d = (Double)Utils.checkNum(decimalsTF.getText(), "display decimals field", null, false, null, new Double(0), true, null, false);
		if(d == null)
			i = 1;
		else
			i = (int)d.doubleValue();

		String s = "0";

		if(i > 0)
			s = s + ".";

		while(i-- > 0)
			s = s + "0";

		unitFmt = new DecimalFormat(s);
	}

	private double avg(final double a, final double b)
	{
		return (a + b) / 2.0;
	}

	private static String makeRow(Object[] o)
	{
		return makeRow(o, false);
	}

	private static String makeRow(Object[] o, boolean header)
	{
		StringBuilder r = new StringBuilder("<tr>");
		String cell = header ? "th" : "td";

		for(int i = 0; i < o.length; i++)
			if(o[i] == null)
				r.append("<" + cell + "></" + cell + ">");
			else
				r.append("<" + cell + ">" + o[i].toString() + "</" + cell + ">");

		r.append("</tr>");
		return r.toString();
	}
}
