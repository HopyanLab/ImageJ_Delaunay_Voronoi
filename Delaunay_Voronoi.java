/*-
 * #%L
 * VIB plugin for Fiji.
 * %%
 * Copyright (C) 2009 - 2024 Fiji developers.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
import delaunay.DelaunayTriangulation;
import delaunay.Pnt;
import delaunay.Simplex;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.gui.StackWindow;
import ij.macro.Interpreter;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.util.Iterator;
import java.util.TreeMap;
import java.lang.Math;

/** Draw Delaunay triangulation or Voronoi Diagram as an overlay. */
public class Delaunay_Voronoi implements PlugIn {

	public final int DELAUNAY = 1;
	public final int VORONOI = 2;
	int mode = DELAUNAY;

	boolean showMeanDistance = false;

	final boolean drawZoom = IJ.getVersion().compareTo("1.37n") >= 0;

	@Override
	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		if (imp == null)
			return;
		GenericDialog gd = new GenericDialog("Delaunay/Voronoi parameters");
		gd.addChoice("mode", new String[] { "Delaunay", "Voronoi"},
				"Delaunay");
		gd.addCheckbox("interactive", !Interpreter.isBatchMode());
		gd.addCheckbox("make Delaunay ROI", true);
		gd.addCheckbox("showMeanDistance", false);
		ResultsTable results = Analyzer.getResultsTable();
		gd.addCheckbox("inferSelectionFromParticles",
				imp.getRoi() == null && results != null
				&& results.getColumnIndex("XM")
				!= ResultsTable.COLUMN_NOT_FOUND);
		gd.addCheckbox("export into Results", true);
		gd.addCheckbox("exclude narrow outer triangles", true);
		gd.addNumericField("minimum angle: ", 25);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
	
		mode = gd.getNextChoiceIndex() + 1;
		boolean interactive = gd.getNextBoolean();
		boolean makeROI = gd.getNextBoolean();
		showMeanDistance = gd.getNextBoolean();
		boolean fromParticles = gd.getNextBoolean();
		boolean exportResults = gd.getNextBoolean();
		boolean excludeNarrow = gd.getNextBoolean();
		double minimumAngle = gd.getNextNumber();
		if (fromParticles) {
			Calibration calib = imp.getCalibration();
			int xCol = results.getColumnIndex("XM");
			int yCol = results.getColumnIndex("YM");
			if (xCol < 0 || yCol < 0) {
				IJ.error("You did not select Center of Mass in"
					+ " Analyze>Set Measurements...\n"
					+ "Select it and try again.");
				return;
			}
			float[] x = results.getColumn(xCol);
			if (x == null || x.length == 0) {
				IJ.error("No results found!");
				return;
			}
			float[] y = results.getColumn(yCol);
			int[] xInt = new int[x.length];
			int[] yInt = new int[x.length];
			for (int i = 0; i < x.length; i++) {
				xInt[i] = (int)Math.round((x[i] /
							calib.pixelWidth) -
							calib.xOrigin);
				yInt[i] = (int)Math.round((y[i] /
							calib.pixelHeight) -
							calib.yOrigin);
			}
			imp.setRoi(new PointRoi(xInt, yInt, x.length));
		}

		CustomCanvas cc = new CustomCanvas(imp);
		TreeMap edges  = getEdges(cc.delaunay, cc.inf,
							excludeNarrow, minimumAngle);
		cc.setEdges(edges);
		if (exportResults) {
			if (results == null) {
				results = new ResultsTable();
				Analyzer.setResultsTable(results);
			}
			exportResults(edges, results);
		}

		if (makeROI) {
			imp.setRoi(getRoi(edges));
			imp.updateAndDraw();
			return;
		}

		if (!interactive) {
			cc.drawOverlay(null);
			imp.updateAndDraw();
			return;
		}

		if (imp.getStackSize()>1)
			new StackWindow(imp, cc).addKeyListener(cc);
		else
			new ImageWindow(imp, cc).addKeyListener(cc);
		Roi roi = imp.getRoi();
		if (roi != null)
			// implicitely set the new image canvas
			roi.setImage(imp);
	}

	void exportResults(TreeMap edges, ResultsTable results) {
		if (mode != DELAUNAY) {
			IJ.error("Operation only supported for Delaunay");
			return;
		}
		if (results.getLastColumn() >= 0) {
			if (!IJ.showMessageWithCancel("Clear Results?",
					"May I clear the results table?"))
				return;
			results.reset();
		}
		TreeMap shown = new TreeMap();
		for (Iterator iter = edges.keySet().iterator(); iter.hasNext();) {
			PntPair pair = (PntPair)iter.next();
			addOneResult(shown, pair.a, pair.b, results);
		}
		results.show("Results");
	}

	TreeMap getEdges (DelaunayTriangulation delaunay, double inf,
			boolean excludeNarrow, double minimumAngle){
		// Assemble Map of all edges and outer edges
		TreeMap edges = new TreeMap();
		TreeMap outer = new TreeMap();
		for (Iterator iter = delaunay.iterator();
				iter.hasNext(); ) {
			Simplex triangle = (Simplex)iter.next();
			Iterator iter2 = triangle.iterator();
			Pnt a = (Pnt)iter2.next();
			Pnt b = (Pnt)iter2.next();
			Pnt c = (Pnt)iter2.next();
			if (Math.abs(a.coord(0)) >= inf ||
					Math.abs(b.coord(0)) >= inf ||
					Math.abs(c.coord(0)) >= inf)
				continue;
			PntPair pair = new PntPair(a, b);
			if (!edges.containsKey(pair)) {
				edges.put(pair, null);
				outer.put(pair, null);
			} else {
				outer.remove(pair);
			}
			pair = new PntPair(a, c);
			if (!edges.containsKey(pair)) {
				edges.put(pair, null);
				outer.put(pair, null);
			} else {
				outer.remove(pair);
			}
			pair = new PntPair(b, c);
			if (!edges.containsKey(pair)) {
				edges.put(pair, null);
				outer.put(pair, null);
			} else {
				outer.remove(pair);
			}
		}
		if (excludeNarrow) {
			double maximumCos = Math.cos(Math.PI*minimumAngle/180.0);
			boolean removed = true;
			while(removed) {
				removed = false;
				TreeMap newEdges = (TreeMap)edges.clone();
				TreeMap newOuter = (TreeMap)outer.clone();
				for (Iterator iter = delaunay.iterator();
					iter.hasNext(); ) {
					Simplex triangle = (Simplex)iter.next();
					Iterator iter2 = triangle.iterator();
					Pnt a = (Pnt)iter2.next();
					Pnt b = (Pnt)iter2.next();
					Pnt c = (Pnt)iter2.next();
					if (Math.abs(a.coord(0)) >= inf ||
							Math.abs(b.coord(0)) >= inf ||
							Math.abs(c.coord(0)) >= inf)
						continue;
					PntPair pair1 = new PntPair(a, b);
					PntPair pair2 = new PntPair(b, c);
					PntPair pair3 = new PntPair(a, c);
					if ((!edges.containsKey(pair1)) ||
						(!edges.containsKey(pair2)) ||
						(!edges.containsKey(pair3)))
						continue;
					if ((!outer.containsKey(pair1)) &&
						(!outer.containsKey(pair2)) &&
						(!outer.containsKey(pair3)))
						continue;
					double ux = b.coord(0) - a.coord(0);
					double uy = b.coord(1) - a.coord(1);
					double vx = c.coord(0) - b.coord(0);
					double vy = c.coord(1) - b.coord(1);
					double wx = a.coord(0) - c.coord(0);
					double wy = a.coord(1) - c.coord(1);
					if (outer.containsKey(pair1)) {
						double num1 = -(ux*vx+uy*vy);
						double num2 = -(ux*wx+uy*wy);
						double den1 = (Math.sqrt(Math.pow(ux, 2) + Math.pow(uy, 2)) * (Math.sqrt(Math.pow(vx, 2) + Math.pow(vy, 2))) );
						double den2 = (Math.sqrt(Math.pow(ux, 2) + Math.pow(uy, 2)) * (Math.sqrt(Math.pow(wx, 2) + Math.pow(wy, 2))) );
						double cos1 = num1/den1;
						double cos2 = num2/den2;
						if ((cos1>maximumCos)||(cos2>maximumCos)) {
							removed = true;
							newEdges.remove(pair1);
							newOuter.remove(pair1);
							if ((!outer.containsKey(pair2)) && 
									(edges.containsKey(pair2)))
								newOuter.put(pair2, null);
							if ((!outer.containsKey(pair3)) && 
									(edges.containsKey(pair3)))
								newOuter.put(pair3, null);
						}
					} else
					if (outer.containsKey(pair2)) {
						double num1 = -(ux*vx+uy*vy);
						double num2 = -(vx*wx+vy*wy);
						double den1 = (Math.sqrt(Math.pow(ux, 2) + Math.pow(uy, 2)) * (Math.sqrt(Math.pow(vx, 2) + Math.pow(vy, 2))) );
						double den2 = (Math.sqrt(Math.pow(vx, 2) + Math.pow(vy, 2)) * (Math.sqrt(Math.pow(wx, 2) + Math.pow(wy, 2))) );
						double cos1 = num1/den1;
						double cos2 = num2/den2;
						if ((cos1>maximumCos)||(cos2>maximumCos)) {
							removed = true;
							newEdges.remove(pair2);
							newOuter.remove(pair2);
							if ((!outer.containsKey(pair1)) && 
									(edges.containsKey(pair1)))
								newOuter.put(pair1, null);
							if ((!outer.containsKey(pair3)) && 
									(edges.containsKey(pair3)))
								newOuter.put(pair3, null);
						}
					} else
					if (outer.containsKey(pair3)) {
						double num1 = -(wx*vx+wy*vy);
						double num2 = -(ux*wx+uy*wy);
						double den1 = (Math.sqrt(Math.pow(wx, 2) + Math.pow(wy, 2)) * (Math.sqrt(Math.pow(vx, 2) + Math.pow(vy, 2))) );
						double den2 = (Math.sqrt(Math.pow(ux, 2) + Math.pow(uy, 2)) * (Math.sqrt(Math.pow(wx, 2) + Math.pow(wy, 2))) );
						double cos1 = num1/den1;
						double cos2 = num2/den2;
						if ((cos1>maximumCos)||(cos2>maximumCos)) {
							removed = true;
							newEdges.remove(pair3);
							newOuter.remove(pair3);
							if ((!outer.containsKey(pair2)) && 
									(edges.containsKey(pair2)))
								newOuter.put(pair2, null);
							if ((!outer.containsKey(pair1)) && 
									(edges.containsKey(pair1)))
								newOuter.put(pair1, null);
						}
					}
				}
				edges = newEdges;
				outer = newOuter;
			}
		}
	return edges;
	}

	private static class PntPair implements Comparable {
		Pnt a, b;

		PntPair(Pnt a, Pnt b) {
			if (compare(a, b) > 0) {
				this.a = b;
				this.b = a;
			} else {
				this.a = a;
				this.b = b;
			}
		}

		@Override
		public int compareTo(Object other) {
			PntPair o = (PntPair)other;
			int result = compare(a, o.a);
			if (result == 0)
				result = compare(b, o.b);
			return result;
		}

		public static int compare(Pnt a, Pnt b) {
			double result = a.coord(0) - b.coord(0);
			if (result == 0)
				result = a.coord(1) - b.coord(1);
			return result < 0 ? -1 : result > 0 ? +1 : 0;
		}
	}

	private void addOneResult(TreeMap shown, Pnt a, Pnt b,
			ResultsTable results) {
		PntPair pair = new PntPair(a, b);
		if (shown.containsKey(pair))
			return;
		results.incrementCounter();

		results.addValue("x1", pair.a.coord(0));
		results.addValue("y1", pair.a.coord(1));
		results.addValue("x2", pair.b.coord(0));
		results.addValue("y2", pair.b.coord(1));
		results.addValue("length", Math.sqrt( Math.pow(pair.a.coord(0)-pair.b.coord(0),2) + Math.pow(pair.a.coord(1)-pair.b.coord(1),2) ));
		shown.put(pair, null);
	}

	Roi getRoi(TreeMap edges) {
		if ((edges == null) || (edges.isEmpty()))
			return null;
		if (mode != DELAUNAY) {
			IJ.error("Operation only supported for Delaunay");
			return null;
		}
		int i = 0;
		GeneralPath path = new GeneralPath(GeneralPath.WIND_EVEN_ODD);
		Line2D.Double line = null;
		for (Iterator iter = edges.keySet().iterator(); iter.hasNext();) {
			PntPair pair = (PntPair)iter.next();
			Pnt a = pair.a;
			Pnt b = pair.b;
			int[] x = new int[2];
			int[] y = new int[2];
			x[0] = (int)Math.round(a.coord(0));
			y[0] = (int)Math.round(a.coord(1));
			x[1] = (int)Math.round(b.coord(0));
			y[1] = (int)Math.round(b.coord(1));
			line = new Line2D.Double(x[0], y[0], x[1], y[1]);
			path.append(line, false);
			i++;
		}
		if (i == 0)
			return null;
		return new ShapeRoi(path);
	}

	Roi getRoi(DelaunayTriangulation delaunay, double inf) {
		if (delaunay == null)
			return null;
		if (mode != DELAUNAY) {
			IJ.error("Operation only supported for Delaunay");
			return null;
		}
		int i = 0;
		GeneralPath path = new GeneralPath(GeneralPath.WIND_EVEN_ODD);
		Polygon poly = null;
		for (Iterator iter = delaunay.iterator();
				iter.hasNext(); ) {
			Simplex triangle = (Simplex)iter.next();
			if (mode == DELAUNAY) {
				Iterator iter2 = triangle.iterator();
				Pnt a = (Pnt)iter2.next();
				Pnt b = (Pnt)iter2.next();
				Pnt c = (Pnt)iter2.next();
				if (Math.abs(a.coord(0)) >= inf ||
						Math.abs(b.coord(0)) >= inf ||
						Math.abs(c.coord(0)) >= inf)
					continue;
				int[] x = new int[3];
				int[] y = new int[3];
				x[0] = (int)Math.round(a.coord(0));
				y[0] = (int)Math.round(a.coord(1));
				x[1] = (int)Math.round(b.coord(0));
				y[1] = (int)Math.round(b.coord(1));
				x[2] = (int)Math.round(c.coord(0));
				y[2] = (int)Math.round(c.coord(1));
				poly = new Polygon(x, y, 3);
				path.append(poly, false);
				i++;
			} else {
				return null;
				/*
					TODO:
				Iterator iter2 = delaunay
					.neighbors(triangle).iterator();
				while (iter2.hasNext())
					draw(g, triangle,
							(Simplex)iter2.next());
				*/
			}
		}
		if (i == 0)
			return null;
		if (i == 1)
			return new PolygonRoi(poly, PolygonRoi.POLYGON);
		return new ShapeRoi(path);
	}

	class CustomCanvas extends ImageCanvas implements KeyListener {
		DelaunayTriangulation delaunay;
		final double inf;
		TreeMap edges = new TreeMap();

		CustomCanvas(ImagePlus imp) {
			super(imp);
			inf = imp.getWidth() + imp.getHeight();
			initDelaunay();
			addKeyListener(this);
		}

		@Override
		public void paint(Graphics g) {
			super.paint(g);
			drawOverlay(g);
		}

		void setEdges (TreeMap edges) {
			this.edges = edges;
		}

		void drawOverlay(Graphics g) {
			if (delaunay == null)
				return;
			IJ.error(String.valueOf(this.edges.isEmpty()));
			if (this.edges.isEmpty()) {
				for (Iterator iter = delaunay.iterator();
						iter.hasNext(); ) {
					Simplex triangle = (Simplex)iter.next();
					if (mode == DELAUNAY) {
						Iterator iter2 = triangle.iterator();
						Pnt a = (Pnt)iter2.next();
						Pnt b = (Pnt)iter2.next();
						Pnt c = (Pnt)iter2.next();
						draw(g, a, b);
						draw(g, a, c);
						draw(g, b, c);
					} else {
						Iterator iter2 = delaunay
							.neighbors(triangle).iterator();
						while (iter2.hasNext())
							draw(g, triangle,
								(Simplex)iter2.next());
					}
				}
			} else {
				for (Iterator iter = edges.keySet().iterator();
							iter.hasNext();) {
					PntPair pair = (PntPair)iter.next();
					draw(g, pair.a, pair.b);
				}
			}
		}

		void draw(Graphics g, Pnt a, Pnt b) {
			if (mode != VORONOI && (Math.abs(a.coord(0)) >= inf ||
						Math.abs(b.coord(0)) >= inf))
				return;

			if (g == null) {
				ImageProcessor ip = imp.getProcessor();
				ip.drawLine((int)a.coord(0),
						(int)a.coord(1),
						(int)b.coord(0),
						(int)b.coord(1));
				return;
			}

			double m = magnification;
			double x0 = (a.coord(0) - srcRect.x) * m;
			double y0 = (a.coord(1) - srcRect.y) * m;
			double x1 = (b.coord(0) - srcRect.x) * m;
			double y1 = (b.coord(1) - srcRect.y) * m;
			g.setColor(imp.getRoi().getColor());
			g.drawLine((int)x0, (int)y0, (int)x1, (int)y1);
			if (drawZoom && srcRect.width != imageWidth
					&& g != null) {
				int xOffset = 10, yOffset = 10;
				int w = 64, h = 64;
				if (imageHeight > imageWidth) {
					m = 64.0 / imageHeight;
					w = (int)(imageWidth * m);
				} else {
					m = 64.0 / imageWidth;
					h = (int)(imageHeight * m);
				}
				x0 = a.coord(0) * m + xOffset;
				y0 = a.coord(1) * m + yOffset;
				x1 = b.coord(0) * m + xOffset;
				y1 = b.coord(1) * m + yOffset;
				Shape clip = g.getClip();
				g.setColor(new Color(128, 128, 255));
				g.clipRect(xOffset, yOffset, w, h);
				g.drawLine((int)x0, (int)y0,
						(int)x1, (int)y1);
				g.setClip(clip);
			}
		}

		void draw(Graphics g, Simplex a, Simplex b) {
			draw(g, Pnt.circumcenter((Pnt[])a.toArray(new Pnt[0])),
				Pnt.circumcenter((Pnt[])b.toArray(new Pnt[0])));
		}

		@Override
		public void mouseReleased(MouseEvent e) {
			super.mouseReleased(e);
			initDelaunay();
			repaint();
		}

		@Override
		public void keyTyped(KeyEvent e) {}
		@Override
		public void keyPressed(KeyEvent e) {
			if (e.getKeyCode() == KeyEvent.VK_SPACE) {
				mode = mode == DELAUNAY ? VORONOI : DELAUNAY;
				repaint();
			}
		}
		@Override
		public void keyReleased(KeyEvent e) {}

		public void initDelaunay() {
			delaunay = null;

			Roi roi = imp.getRoi();
			if (roi == null || !(roi instanceof PointRoi))
				return;

			PointRoi r = (PointRoi)roi;
			Rectangle rect = r.getBounds();
			int n = r.getNCoordinates();
			int[] x = r.getXCoordinates();
			int[] y = r.getYCoordinates();

			Simplex initial = new Simplex(new Pnt[] {
					new Pnt(-inf, -inf),
					new Pnt(-inf, 5 * inf),
					new Pnt(5 * inf, -inf)});
			delaunay = new DelaunayTriangulation(initial);
			for (int i = 0; i < n; i++)
				delaunay.delaunayPlace(new Pnt(x[i] + rect.x,
							y[i] + rect.y));

			if (showMeanDistance && mode == DELAUNAY)
				showMeanAndVariance();
		}

		double pixelWidth, pixelHeight;
		double mean, variance;
		int total;
		
		private void addToMean(Pnt a, Pnt b) {
			if (Math.abs(a.coord(0)) >= inf ||
					Math.abs(b.coord(0)) >= inf)
				return;
			double x = (b.coord(0) - a.coord(0)) * pixelWidth;
			double y = (b.coord(1) - a.coord(1)) * pixelHeight;
			double d2 = x * x + y * y;
			mean += Math.sqrt(d2);
			variance += d2;
			total++;
		}

		public void showMeanAndVariance() {
			Calibration calib = imp.getCalibration();
			pixelWidth = calib.pixelWidth;
			pixelHeight = calib.pixelHeight;

			mean = variance = total = 0;

			for (Iterator iter = delaunay.iterator();
					iter.hasNext(); ) {
				Simplex triangle = (Simplex)iter.next();
				Iterator iter2 = triangle.iterator();
				Pnt a = (Pnt)iter2.next();
				Pnt b = (Pnt)iter2.next();
				Pnt c = (Pnt)iter2.next();
				addToMean(a, b);
				addToMean(b, c);
				addToMean(c, a);
			}

			if (total > 0) {
				mean /= total;
				variance /= total;
				variance -= mean * mean;

				ResultsTable rt = Analyzer.getResultsTable();
				if (rt == null) {
					rt = new ResultsTable();
					Analyzer.setResultsTable(rt);
				}

				rt.incrementCounter();
				rt.addValue("Mean Distance", mean);
				rt.addValue("Variance", variance);
				rt.show("Results");
			}
		}
	}
}
