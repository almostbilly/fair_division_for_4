import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.StandardCategoryToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.CombinedDomainCategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

import java.awt.*;

public class BarChart {
    private static final long serialVersionUID = 1L;
    public CategoryDataset dataset1, dataset2;
    public BarChart() { }

    public void setDataset1(CategoryDataset dataset1) {
        this.dataset1 = dataset1;
    }

    public void setDataset2(CategoryDataset dataset2) {
        this.dataset2 = dataset2;
    }

    public JFreeChart createChart(CategoryDataset dataset1, CategoryDataset dataset2) {
        NumberAxis rangeAxis1 = new NumberAxis("Стоимость после 1го деления");
        NumberAxis rangeAxis2 = new NumberAxis("Стоимость в конце");
        rangeAxis1.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        rangeAxis2.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        BarRenderer renderer = new BarRenderer();
        renderer.setBaseToolTipGenerator(new StandardCategoryToolTipGenerator());
        CategoryPlot subplot1 = new CategoryPlot(dataset1, null, rangeAxis1, renderer);
        CategoryPlot subplot2 = new CategoryPlot(dataset2, null, rangeAxis2, renderer);
        CategoryAxis domainAxis = new CategoryAxis("Оценка");
        CombinedDomainCategoryPlot plot = new CombinedDomainCategoryPlot(domainAxis);
        plot.add(subplot1);
        plot.add(subplot2);
        JFreeChart chart = new JFreeChart("Оценка агентов", plot);
        return chart;
    }

    public void showChart() {
        JFreeChart chart = createChart(dataset1, dataset2);
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(1024,768));
        ApplicationFrame frame = new ApplicationFrame("Plots");
        frame.setContentPane(chartPanel);
        frame.pack();
        RefineryUtilities.centerFrameOnScreen(frame);
        frame.setVisible(true);
    }

}
