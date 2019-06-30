using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Forms.DataVisualization.Charting;
using System.Drawing;
using System.Diagnostics;


namespace SignalProcessingTool
{
    using System.Collections.Generic;
    using OxyPlot;
    using System.Diagnostics;
    using OxyPlot.Series;
    using OxyPlot.Axes;
    using System.Collections.ObjectModel;

    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();

        }

        private void BtnCalcStart_Click(object sender, RoutedEventArgs e)
        {
            AcousticModel model1 = new AcousticModel();
            
            CalculateImpulse();
        }

        public void CalculateImpulse()
        {

            
        }

    }

    public class MainViewModel
    {
        public MainViewModel()
        {
            this.Points = new Collection<PointClass>();
            const int N = 4096;

            for (double i = 0; i < 1; i = i + 0.001)
            {               
                this.Points.Add(
                new PointClass
                {
                    xPoint = i,
                    yPoint = Math.Sin(i*10),
                    yPointMax = Math.Sin(i * 10)/2
                    
                });
            }
        }

        public Collection<PointClass> Points { get; private set; }

        public string Subtitle { get; set; }

    }

    public class AcousticModel
    {
        public double Func { get; set; }
        public double Time { get; set; }
        public double Mass { get; set; }
        public double Speed { get; set; }
        public double Density { get; set; }
        public double Area { get; set; }
        public double PipeLength { get; set; }
        public double Tau { get; set; }
        public double xCoordinate { get; set; }
        public double modNumber { get; set; }
        public double modFrequency { get; set; }
        public double Delta { get; set; }
        public double A1 { get; set; }
        public double A2 { get; set; }
        public double A3 { get; set; }
    }

    public class PointClass
    {
        public double xPoint { get; set; }
        public double yPoint { get; set; }
        public double yPointMax { get; set; }
    }
}
