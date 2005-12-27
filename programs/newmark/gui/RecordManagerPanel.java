/*
 * RecordManagerPanel.java - manage records
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

/* $Id$ */

package newmark.gui;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.border.*;
import java.io.File;
import java.util.Vector;
import org.jfree.data.xy.*;
import org.jfree.chart.*;
import org.jfree.chart.axis.*;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.Range;
import newmark.*;
import newmark.analysis.*;

class RecordManagerPanel extends JPanel implements ActionListener
{
	NewmarkTabbedPane parent;

	NewmarkTable table;

	JComboBox eqList = new JComboBox();
	JButton graph = new JButton("Graph");
	JButton delete = new JButton("Delete selected record(s) from database");

	JButton save = new JButton("Save changes");

	JTextField modFile = new JTextField();
	JTextField modEq = new JTextField();
	JTextField modRec = new JTextField();
	JTextField modDI = new JTextField();
	JTextField modLoc = new JTextField();
	JTextField modMag = new JTextField();
	JTextField modOwn = new JTextField();
	JTextField modEpi = new JTextField();
	JTextField modLat = new JTextField();
	JTextField modFoc = new JTextField();
	JTextField modLng = new JTextField();
	JTextField modRup = new JTextField();
	JComboBox  modSite = new JComboBox(NewmarkTable.SiteClassArray);
	JComboBox  modMech = new JComboBox(NewmarkTable.FocMechArray);

	ButtonGroup TypeGroup = new ButtonGroup();
	JRadioButton typeTime = new JRadioButton("Time Series", true);
	JRadioButton typeFourier = new JRadioButton("Fourier Amplitude Spectrum");
	JRadioButton typeSpectra = new JRadioButton("Response Spectra");

	JComboBox spectraCB = new JComboBox();

	JTextField spectraDamp = new JTextField("0", 5);
	JTextField spectraIncr = new JTextField("0.01", 5);
	JTextField spectraHigh = new JTextField("15.0", 5);

	final public static String spectraAccStr = "Absolute-Acceleration";
	final public static String spectraVelStr = "Relative-Velocity";
	final public static String spectraDisStr = "Relative-Displacement";

	JButton add = new JButton("Add record(s)...");

	JTabbedPane managerTP = new JTabbedPane();

	public RecordManagerPanel(NewmarkTabbedPane parent) throws Exception
	{
		this.parent = parent;

		table = new NewmarkTable(false);

		modFile.setEditable(false);

		ListSelectionModel recordSelect = table.getSelectionModel();
		recordSelect.addListSelectionListener(new ListSelectionListener()
		{
			public void valueChanged(ListSelectionEvent event)
			{
				try
				{
					recordSelect(event);
				}
				catch (Exception e)
				{
					Utils.catchException(e);
				}
			}
		});

		TypeGroup.add(typeTime);
		TypeGroup.add(typeFourier);
		TypeGroup.add(typeSpectra);

		spectraCB.addItem(spectraAccStr);
		spectraCB.addItem(spectraVelStr);
		spectraCB.addItem(spectraDisStr);

		Utils.addEQList(eqList, Boolean.TRUE);
		eqList.setActionCommand("eqListChange");
		eqList.addActionListener(this);

		save.setActionCommand("save");
		save.addActionListener(this);

		delete.setActionCommand("delete");
		delete.addActionListener(this);

		graph.setActionCommand("graph");
		graph.addActionListener(this);

		managerTP.addTab("Modify Record", createModifyPanel());
		managerTP.addTab("Graphing Options", createGraphPanel());

		setLayout(new BorderLayout());
		add(BorderLayout.NORTH, createNorthPanel());
		add(BorderLayout.CENTER, table);
		add(BorderLayout.SOUTH, managerTP);

		recordClear();
	}

	private JPanel createNorthPanel()
	{
		JPanel panel = new JPanel(new BorderLayout());

		Vector list = new Vector();
		list.add(new JLabel("Display records from: "));
		list.add(eqList);
		panel.add(BorderLayout.WEST, GUIUtils.makeRecursiveLayoutRight(list));

		list = new Vector();
		list.add(graph);
		list.add(delete);
		panel.add(BorderLayout.EAST, GUIUtils.makeRecursiveLayoutRight(list));

		return panel;
	}

	private JPanel createModifyPanel()
	{
		JPanel panel = new JPanel(new BorderLayout());

		JPanel north = new JPanel(new BorderLayout());

		JLabel label = new JLabel("Modify record:");
		label.setFont(GUIUtils.headerFont);

		JPanel file = new JPanel(new BorderLayout());
		file.add(BorderLayout.WEST, new JLabel("File Location "));
		file.add(BorderLayout.CENTER, modFile);

		north.add(BorderLayout.WEST, label);
		north.add(BorderLayout.EAST, save);
		north.add(BorderLayout.SOUTH, file);

		panel.add(BorderLayout.NORTH, north);

		JPanel south = new JPanel(new GridLayout(0, 4));
		south.add(new JLabel("Earthquake name"));
		south.add(modEq);
		south.add(new JLabel("     Record name"));
		south.add(modRec);
		south.add(new JLabel("Digitization Interval (s)"));
		south.add(modDI);
		south.add(new JLabel("     Location [optional]"));
		south.add(modLoc);
		south.add(new JLabel("Moment Magnitude [optional]"));
		south.add(modMag);
		south.add(new JLabel("     Station owner [optional]"));
		south.add(modOwn);
		south.add(new JLabel("Epicentral distance (km) [optional]"));
		south.add(modEpi);
		south.add(new JLabel("     Latitude [optional]"));
		south.add(modLat);
		south.add(new JLabel("Focal distance (km) [optional]"));
		south.add(modFoc);
		south.add(new JLabel("     Longitude [optional]"));
		south.add(modLng);
		south.add(new JLabel("Rupture distance (km) [optional]"));
		south.add(modRup);
		south.add(new JLabel("     Site Class [optional]"));
		south.add(modSite);
		south.add(new JLabel("Focal Mechanism [optional]"));
		south.add(modMech);

		panel.add(BorderLayout.SOUTH, south);

		return panel;
	}

	public JPanel createGraphPanel()
	{
		JPanel panel = new JPanel();

		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		JLabel label;

		panel.setLayout(gridbag);

		int x = 0;
		int y = 0;

		c.anchor = GridBagConstraints.NORTHWEST;

		c.gridx = x++;
		c.gridy = y++;
		gridbag.setConstraints(typeTime, c);
		panel.add(typeTime);

		c.gridy = y++;
		gridbag.setConstraints(typeFourier, c);
		panel.add(typeFourier);

		c.gridy = y++;
		gridbag.setConstraints(typeSpectra, c);
		panel.add(typeSpectra);

		c.gridx = x++;
		c.gridy = y++;
		gridbag.setConstraints(spectraCB, c);
		panel.add(spectraCB);

		x = 1;
		c.gridx = x++;
		c.gridy = y++;
		label = new JLabel("Damping (%)");
		gridbag.setConstraints(label, c);
		panel.add(label);

		c.gridx = x++;
		label = new JLabel("Frequency Intrement (Hz)");
		gridbag.setConstraints(label, c);
		panel.add(label);

		c.gridx = x++;
		label = new JLabel("High Frequency (Hz)");
		gridbag.setConstraints(label, c);
		panel.add(label);

		x = 1;
		c.gridx = x++;
		c.gridy = y++;
		gridbag.setConstraints(spectraDamp, c);
		panel.add(spectraDamp);

		c.gridx = x++;
		gridbag.setConstraints(spectraIncr, c);
		panel.add(spectraIncr);

		c.gridx = x++;
		gridbag.setConstraints(spectraHigh, c);
		panel.add(spectraHigh);

		return panel;
	}

	public void actionPerformed(java.awt.event.ActionEvent e)
	{
		try
		{
			String command = e.getActionCommand();
			System.out.println(command);
			if(command.equals("delete"))
			{
				int n = JOptionPane.showConfirmDialog(this, "Do you want to delete these records?", "Delete?", JOptionPane.YES_NO_OPTION);
				if(n != JOptionPane.YES_OPTION)
					return;

				table.deleteSelected(true);
			}
			else if(command.equals("eqListChange"))
			{
				if(Utils.locked())
					return;

				boolean isEq = true;
				int index = eqList.getSelectedIndex();
				for(int i = 0; i < eqList.getItemCount(); i++)
				{
					if(((String)eqList.getItemAt(i)).equals(" -- Groups -- "))
					{
						if(i == index)
							return;
						else if(index > i)
							isEq = false;
						break;
					}
				}

				Utils.getDB().runUpdate("update data set select1=0 where select1=1");

				if(isEq)
				{
					String where;

					if(index == 0) // "All earthquakes"
						where = "";
					else
						where = "where eq='" + (String)eqList.getSelectedItem() + "'";

					Utils.getDB().runUpdate("update data set select1=1 " + where);
				}
				else
				{
					Object[][] res = Utils.getDB().runQuery("select record,analyze from grp where name='" + (String)eqList.getSelectedItem() + "'");

					for(int i = 1; i < res.length; i++)
						Utils.getDB().runUpdate("update data set select1=1 where id=" + res[i][0].toString());
				}

				table.setModel(NewmarkTable.REFRESH);
				recordClear();
			}
			else if(command.equals("graph"))
			{
				int row = table.getSelectedRow();
				if(row == -1)
					return;

				String eq = table.getModel().getValueAt(row, 0).toString();
				String record = table.getModel().getValueAt(row, 1).toString();

				Object[][] res = null;
				res = Utils.getDB().runQuery("select path,digi_int from data where eq='" + eq + "' and record='" + record + "'");
				if(res == null) return;

				String path = res[1][0].toString();
				double di = Double.parseDouble(res[1][1].toString());

				File f = new File(path);
				if(f.canRead() == false)
				{
					GUIUtils.popupError("Cannot read or open the file " + path + ".");
					return;
				}

				DoubleList dat = new DoubleList(path);

				if(dat.bad())
				{
					GUIUtils.popupError("Invalid data at data point " + dat.badEntry() + " in " + path + ".");
					return;
				}

				XYSeries xys = new XYSeries("");
				dat.reset();

				String xAxis = null, yAxis = null, title = "";

				if(typeTime.isSelected()) command = "graphTime";
				else if(typeFourier.isSelected()) command = "graphFourier";
				else if(typeSpectra.isSelected()) command = "graphSpectra";

				if(command.equals("graphTime"))
				{
					title = "Time Series";
					xAxis = "Time (s)";
					yAxis = "Acceleration (cm/s/s)";
					Double val;
					double last1 = 0, last2 = 0, current;
					double diff1, diff2;
					double time = 0, timeStor = 0, td;
					int perSec = 50;
					double interval = 1.0 / (double)perSec;

					// add the first point
					if((val = dat.each()) != null)
					{
						xys.add(time, val);
						time += di;
						last2 = val.doubleValue();
					}

					// don't add the second point, but update the data
					if((val = dat.each()) != null)
					{
						time += di;
						last1 = val.doubleValue();
					}

					while((val = dat.each()) != null)
					{
						td = time - di;
						current = val.doubleValue();

						diff1 = last1 - current;
						diff2 = last1 - last2;

						if(
							(diff1 <= 0 && diff2 <= 0) ||
							(diff1 >= 0 && diff2 >= 0) ||
							(td >= (timeStor + interval)))
						{
							xys.add(td, last1);
							timeStor = td;
						}

						last2 = last1;
						last1 = current;
						time += di;
					}
				}
				else if(command.equals("graphFourier"))
				{
					title = "Fourier Amplitude Spectrum";
					xAxis = "Frequency (Hz)";
					yAxis = "Fourier Amplitude (cm/sec)";
					double[] arr = new double[dat.size()];

					Double temp;
					for(int i = 0; (temp = dat.each()) != null; i++)
						arr[i] = temp.doubleValue();

					double[][] fft = ImportRecords.fftWrap(arr, di);

					double df = 1.0 / ((double)(arr.length) * di);

					double current;
					double freq = 0;
					int step, i;

					for(i = 0; i < arr.length;)
					{
						// don't graph anything below 0.1 hz
						if(freq > 0.1)
							xys.add(freq, Math.sqrt(Math.pow(fft[i][0], 2) + Math.pow(fft[i][1], 2)));

						// throw out 3/4 of the points above 10 hz
						if(freq > 10.0)
							step = 4;
						else
							step = 1;

						i += step;
						freq = i * df;
					}
				}
				else if(command.equals("graphSpectra"))
				{
					int index = spectraCB.getSelectedIndex();

					title = spectraCB.getSelectedItem().toString() + " Response Spectra";
					xAxis = "Period (s)";
					yAxis = "Response";
					double[] arr = new double[dat.size()];

					Double temp;
					for(int i = 0; (temp = dat.each()) != null; i++)
						arr[i] = temp.doubleValue();

					double periodMax = Double.parseDouble(spectraHigh.getText());
					double interval = Double.parseDouble(spectraIncr.getText());
					double damp = Double.parseDouble(spectraDamp.getText()) / 100.0;

					double[] z;

					// don't start with 0, LogarithmicAxis doesn't allow it
					for(double p = interval; p < periodMax; p += interval)
					{
						z = ImportRecords.cmpmax(arr, 2.0 * Math.PI / p, damp, di);
						xys.add(p, z[index]);
					}
				}

				title += ": " + eq + " - " + record;
				XYSeriesCollection xysc = new XYSeriesCollection(xys);

				JFreeChart chart = ChartFactory.createXYLineChart(title, xAxis, yAxis, xysc, org.jfree.chart.plot.PlotOrientation.VERTICAL, false, true, false);

				if(command.equals("graphFourier") || command.equals("graphSpectra"))
				{
					chart.getXYPlot().setRangeAxis(new LogarithmicAxis(yAxis));
					chart.getXYPlot().setDomainAxis(new LogarithmicAxis(xAxis));
				}

				ChartFrame frame = new ChartFrame(title, chart);
				frame.pack();
				frame.setLocationRelativeTo(null);
				frame.setVisible(true);
			}
			else if(command.equals("save"))
			{
				int n = JOptionPane.showConfirmDialog(this,"Are you sure you want to modify these records?", "Are you sure?", JOptionPane.YES_NO_OPTION);
				if(n != JOptionPane.YES_OPTION)
					return;

				String error = AddRecordsPanel.manipRecord(false,
					modFile.getText(),
					modEq.getText(),
					modRec.getText(),
					modDI.getText(),
					Utils.nullify(modMag.getText()),
					Utils.nullify(modEpi.getText()),
					Utils.nullify(modFoc.getText()),
					Utils.nullify(modRup.getText()),
					modMech.getSelectedItem().toString(),
					modLoc.getText(),
					modOwn.getText(),
					Utils.nullify(modLat.getText()),
					Utils.nullify(modLng.getText()),
					modSite.getSelectedItem().toString()
				);

				if(!error.equals(""))
					GUIUtils.popupError(error);

				table.setModel(NewmarkTable.REFRESH);
			}
		}
		catch (Exception ex)
		{
			Utils.catchException(ex);
		}
	}

	public void recordSelect(ListSelectionEvent e)
	{
		if(e == null) return;

		if (e.getValueIsAdjusting()) return;

		ListSelectionModel lsm = (ListSelectionModel)e.getSource();

		if (!lsm.isSelectionEmpty())
		{
			int selectedRow = lsm.getMinSelectionIndex();

			String eq = table.getModel().getValueAt(selectedRow, 0).toString();
			String record = table.getModel().getValueAt(selectedRow, 1).toString();

			Object[][] res = null;

			try
			{
				res = Utils.getDB().runQuery("select path, eq, record, digi_int, location, mom_mag, owner, epi_dist, latitude, foc_dist, longitude, rup_dist, class, foc_mech, change from data where eq='" + eq + "' and record='" + record + "'");
			}
			catch(Exception ex)
			{
				Utils.catchException(ex);
			}
			boolean enable = false;
			if(res != null)
			{
				if(res.length > 1)
				{
					int incr = 0;
					modFile.setText(res[1][incr++].toString());
					modEq.setText(res[1][incr++].toString());
					modRec.setText(res[1][incr++].toString());
					modDI.setText(res[1][incr++].toString());
					modLoc.setText(res[1][incr++].toString());
					modMag.setText(Utils.shorten(res[1][incr++]));
					modOwn.setText(res[1][incr++].toString());
					modEpi.setText(Utils.shorten(res[1][incr++]));
					modLat.setText(Utils.shorten(res[1][incr++]));
					modFoc.setText(Utils.shorten(res[1][incr++]));
					modLng.setText(Utils.shorten(res[1][incr++]));
					modRup.setText(Utils.shorten(res[1][incr++]));
					modSite.setSelectedItem(NewmarkTable.SiteClassArray[Integer.parseInt(Utils.shorten(res[1][incr++].toString()))]);
					modMech.setSelectedItem(NewmarkTable.FocMechArray[Integer.parseInt(Utils.shorten(res[1][incr++].toString()))]);
					enable = res[1][incr++].toString().equals("1") ? true : false;
				}
				else
				{
					recordClear();
				}
			}
			else
			{
				recordClear();
			}
			recordEnable(enable);
		}
	}

	public void recordEnable(boolean b)
	{
		modEq.setEditable(b);
		modRec.setEditable(b);
		modDI.setEditable(b);
		modLoc.setEditable(b);
		modMag.setEditable(b);
		modOwn.setEditable(b);
		modEpi.setEditable(b);
		modLat.setEditable(b);
		modFoc.setEditable(b);
		modLng.setEditable(b);
		modRup.setEditable(b);
		modSite.setEnabled(b);
		modMech.setEnabled(b);
		save.setEnabled(b);
	}

	public void recordClear()
	{
		modFile.setText("");
		modEq.setText("");
		modRec.setText("");
		modDI.setText("");
		modLoc.setText("");
		modMag.setText("");
		modOwn.setText("");
		modEpi.setText("");
		modLat.setText("");
		modFoc.setText("");
		modLng.setText("");
		modRup.setText("");
		modSite.setSelectedItem("");
		modMech.setSelectedItem("");

		recordEnable(false);
	}
}
