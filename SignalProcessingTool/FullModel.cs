using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SignalProcessingTool
{
    /// <summary>
    /// A full model of acoustic signal
    /// </summary>
    public class FullAcousticModel
    {        
        public double Diameter { get; set; }
        public double Thickness { get; set; }
        public double Density { get; set; }
        public double PipeLength { get; set; }
        public double XCoordinate { get; set; }
        public double A1 { get; set; }
        public double A2 { get; set; }
        public double A3 { get; set; }
        public double SignalDuration { get; set; }
        public double Tc { get; set; }
        public double Ti { get; set; }
        public double Fd { get; set; }
        public double ModNumber { get; set; }
        public double[] Cn { get; set; }
        public double[] Delta { get; set; }
        public double[] W { get; set; }
        public double[] Oi { get; set; }
        public double[] Sum { get; set; }
        public double E { get; set; }
        public double J { get; set; }
        public double Mass { get; set; }
        public double A4 { get; set; }
        public double Pi { get; set; }

        public virtual void ComputeModel()
        {
            int k = 0;

            for (int i = 1; i < ModNumber; i = i + 2)
            {
                Cn[k] = Pi / 2 * (2 * i + 1);
                Delta[k] = Math.Sqrt(0.5 * (A1 * Math.Pow(i, 4) * Math.Pow(Pi, 4) / Math.Pow(PipeLength, 4) + A2));
                W[k] = Math.Sqrt(A3 + A4 * Math.Pow(i, 4) * Math.Pow(Pi, 4) / Math.Pow(PipeLength, 4));
                Oi[k] = Math.Sqrt(Math.Pow(W[k], 2) - Math.Pow(Delta[k], 2));
                k++;
            }

            k = 0;

            double m1, m2, m3, rotator;

            for (int i = 1; i < SignalDuration * Fd; i++)
            {
                for (int j = 1; j < ModNumber; j = j + 2)
                {
                    rotator = Math.Pow(-1, (j - 1) / 2);

                    m1 = 2 * Math.Exp(Delta[k] * Tc) * Delta[k] * Pi / Tc;

                    m2 = Pi * Math.Exp(Delta[k] * Tc) / Tc / Oi[k] * (2 * Math.Pow(Delta[k], 2) - Math.Pow(W[k], 2) + Math.Pow(Pi, 2) / Math.Pow(Tc, 2));

                    m3 = Oi[k] * Tc;

                    Sum[i] = Sum[i] + rotator * Math.Exp(-Delta[k] * Ti) * Math.Sin(j * Pi * XCoordinate / PipeLength) /
                    (Math.Pow((Math.Pow(W[k], 2) - Math.Pow(Pi, 2) / Math.Pow(Tc, 2)), 2) + 4 * Math.Pow(Pi, 2) * Math.Pow(Delta[k], 2) / Math.Pow(Tc, 2)) *
                    ((m1 * Math.Cos(m3) - m2 * Math.Sin(m3) + 2 * Pi * Delta[k] / Tc) * Math.Cos(Oi[k] * Ti) + (m2 * Math.Cos(m3) + m1 * Math.Sin(m3) +
                    Pi / Tc / Oi[k] * (2 * Math.Pow(Delta[k], 2) - Math.Pow(W[k], 2) + Math.Pow(Pi, 2) / Math.Pow(Tc, 2))) * Math.Sin(Oi[k] * Ti));

                    k++;
                }

                k = 0;

                Ti = Ti + (double)(1 / Fd);
            }
        }
    }
}
