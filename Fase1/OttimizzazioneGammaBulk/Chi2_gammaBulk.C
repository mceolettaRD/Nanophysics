{
	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	const double e0 = 9.84; //epsilon0
	const double hbar = 6.582119569e-16; //eV/s
	//tutte le energie sono in elettronvolt
	const double Ep = 9.010; //frequenza plasma bulk
	const double Ec = 2.4; //parametri per la correzione logistica
	const double Delta = 0.17;
	const double A = 5.6; //adimensionale
	const double lmin = 700; //intervallo in cui generare la dielettrica di bulk, in nm
	const double lmax = 900;
	const double stepping = 1; //step per la generazione, in nm

	double n = 0; //counter per grafico
	double gammab = 0; //gammabulk materiale, verrà variato qui solo inizzializzazione failsafe

	struct diel {
		double l, E, e1, e2;
	};

	//lettura e storage di JC
	vector<diel> dataJC;
	ifstream inp;
	string fil="epsiJC_cut.txt"; //già tagliata da 700 a 900
   	inp.open(fil);

	double l,e1,e2,E;
   	while(inp >> l >> E >> e1 >> e2)
   	{
   		diel d;
   		d.l = l;
   		d.E = E;
   		d.e1 = e1;
   		d.e2 = e2;
   		dataJC.push_back(d);
   	}
   	cout << "dataJCSize = " << dataJC.size() << endl;

   	TGraph *chigd = new TGraph();

   	for (double b = 0.061; b <= 0.07; b+= 0.0001) //ciclo sui valori di gammaBulk
   	{
   		gammab = b;
   		cout << "valore gammab = " << gammab << endl;

   		vector<diel> dataBulk; //generazione del gammabulk

		for (double i = lmin; i <= lmax; i+=stepping)
		{
			diel d;

			d.l = i;
			d.E = 1239.84193/(i);

			d.e1 = e0 - pow(Ep,2)/(pow(d.E, 2)+pow(gammab, 2));
			d.e2 = pow(Ep,2)*gammab/(d.E*(pow(d.E,2)+pow(gammab,2))) + A/(1+exp(-1*(d.E-Ec)/Delta));

			dataBulk.push_back(d);
		}

		//calcolo del chi2

		double chi2 = 0;

		for (int i = 0; i < min(dataBulk.size(), dataJC.size()); ++i)
		{
			chi2 += pow(dataJC[i].e2 - dataBulk[i].e2,2)/dataJC[i].e2;
		}
		cout << "valore chi2 = " << chi2 << endl;
		chigd->SetPoint(n++,gammab,chi2);
	}

	TF1 *parab = new TF1("para","pol2");
	chigd->Fit(parab,"M");

	double vertice = -1*parab->GetParameter(1)/(2*parab->GetParameter(2));
	double a = parab->GetParameter(2);
	double b = parab->GetParameter(1);
	double c = parab->GetParameter(0);
	double p = vertice+1;
	double clLo = (-1*b - sqrt(b*b -4*a*(c-p)))/(2*a);
	double clHi = (-1*b + sqrt(b*b -4*a*(c-p)))/(2*a);
	
	cout << "il vertice = " << vertice << " CL" << clLo << ":" << clHi << endl;

	chigd->SetMarkerStyle(3);
	chigd->GetXaxis()->SetTitle("#Gamma_bulk");
	chigd->GetYaxis()->SetTitle("#chi^2");
	chigd->Draw("ap");
	c1->Update();

	TLine *llo = new TLine(clLo,c1->GetUymin(),clLo,c1->GetUymax());
	TLine *lhi = new TLine(clHi,c1->GetUymin(),clHi,c1->GetUymax());
	TLine *lm = new TLine(vertice,c1->GetUymin(),vertice,c1->GetUymax());
	llo->Draw();
	lhi->Draw();
	c1->SetLineColor(kRed);
	lm->Draw();
}