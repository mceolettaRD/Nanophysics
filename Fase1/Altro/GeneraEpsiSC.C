{
	const double e0 = 9.84; //epsilon0
	const double hbar = 6.582119569e-16; //eV/s

	//tutte le energie sono in elettronvolt
	const double gammab = 0.072; //gammabulk materiale
	const double Ep = 9.010; //energia plasma bulk

	const double Ec = 2.4; //parametri per la correzione logistica
	const double Delta = 0.17;
	const double A = 5.6; //adimensionale, riscalamento parte logistica

	/* EFFETTI DI SCALA, CALCOLO GAMMA(r) */
	const double r = 10; //raggio particelle in nm (monodisperso...)
	const double a = TMath::Pi()/4; //costante di riscalamento della legge di scala per il fattore gamma (~1)
	const double vf = 1.40e15; //velocità di Fermi degli elettroni di bulk 1.40e15 nm/s, divido per nm
	const double gammar = gammab + a*vf/r*hbar;
	cout << gammar << endl;

	const double lmin = 200; //intervallo in cui generare la dielettrica, in nm, 467nm èil taglio per il valore costante...
	const double lmax = 900;
	const double stepping = 1; //step per la generazione, in nm

	struct diel {
		double l, e1, e2;
	};

	vector<diel> dataset;

	for (double i = lmin; i <= lmax; i+=stepping)
	{
		diel d;

		d.l = i;
		double E = 1239.84193/(i);

		d.e1 = e0 - pow(Ep,2)/(pow(E, 2)+pow(gammar, 2));
		d.e2 = pow(Ep,2)*gammar/(E*(pow(E,2)+pow(gammar,2))) + A/(1+exp(-1*(E-Ec)/Delta));
		printf("%.1f,%f : %f:%f\n",d.l, E,d.e1,d.e2);
		dataset.push_back(d);
	}

	ofstream outf("epsiSC.txt");

   	for(int i = 0; i < dataset.size(); i++)
   	{

   		outf << dataset[i].l << " " << 1239.84193/dataset[i].l << " " << dataset[i].e1 << " "<< dataset[i].e2 << endl;
   	}

   	outf.close();
}