{
	//const double epsim = 2.182; //fisso valore di epsim a quello di frolich
	const double e0 = 9.84;
	const double gammab = 0.0654926;
	const double Ep = 9.010;
	const double vf = 1.40e15; //in nm/s
	const double hbar = 6.582119569e-16; //eV/s
	const double epsim = 2;

	struct diel {
		double l, E, e1, e2;
	};

	struct sigma {
	double l, sig = 0;
	};

	struct chi {
	double a1, ecc;
	double chi2 = 0;
	};

	vector<diel> dataJC;
	ifstream inp1;
	string fil1="epsiJC.txt";
	inp1.open(fil1);
	double l,e1,e2,E;
	while(inp1 >> l >> E >> e1 >> e2)
	{
		diel d;
		d.l = l;
		d.E = E;
		d.e1 = e1;
		d.e2 = e2;
		dataJC.push_back(d);
	}
	inp1.close();


	vector<sigma> assSpet;
	ifstream inp2;
	string fil2="AuNP_rev.txt";
	inp2.open(fil2);
	double sig;
	while(inp2 >> l >> sig)
	{
		sigma d;
		d.l = l;
		d.sig = sig;
		assSpet.push_back(d);
	}
	inp2.close();

	//calcolo del valore medio all'inizio dello spettro
	double meanSpet = 0;
	int nSpet = 0;
	for (int i = 0; i < assSpet.size(); ++i)
	{
		if(assSpet[i].l >= 400 && assSpet[i].l <= 420)
		{
			//printf("assSpet = %f,%f \n",assSpet[i].l, assSpet[i].sig);
			meanSpet += assSpet[i].sig;
			nSpet++;
		}
	}
	meanSpet /= nSpet;

	//teoria di Gans, correzione usando semiasse maggiore
	vector<chi> chiTable;
	for(double a1 = 0.5; a1 <= 20; a1 += 0.01)
	{
		 cout << a1 << endl;
		 for (double ecc = 0.; ecc <= 1; ecc += 0.001)
		 {
		 	double a23 = a1*sqrt(1-ecc*ecc);
			double vol = 4*TMath::Pi()*(a1*a23*a23)/3; //volume NP
			//printf("%f, %f \n",r,epsim);

			double L1 = (1-pow(ecc,2))/(pow(ecc, 2))*(log((1+ecc)/(1-ecc))/(2*ecc)-1);
			double L23 = (1-L1)/2;
			double L[] = {L1, L23, L23};
			

			chi entry; // costruisco l'entrata del chiTable e assegno valori dell'iterazione
			entry.a1 = a1;
			entry.ecc = ecc;

			double gammar = gammab + vf/(4*a1)*hbar;

			vector<diel> dataJC_cor;
			for(int i = 0; i < dataJC.size(); i++)
			{
				diel d;
				d.l = dataJC[i].l;
				d.E = dataJC[i].E;
				double beta = (pow(Ep, 2)/d.E)*1/(pow(d.E, 2)-pow(gammar, 2));
				d.e1 = dataJC[i].e1 + pow(Ep,2)*(1/(pow(E,2)+pow(gammab,2)) - 1/(pow(E,2)+pow(gammar,2))); //- beta*d.E
				d.e2 = dataJC[i].e2 - (pow(Ep,2)/E)*( gammab/(pow(E,2)+pow(gammab,2)) - gammar/(pow(E,2)+pow(gammar,2)) ); //+ beta*gammar
				dataJC_cor.push_back(d);
			}

			//calcolo la sigma di JC
			vector<sigma> extJC;
			for(int i = 0; i < dataJC_cor.size(); i++)
			{
				sigma s;
				s.l = dataJC_cor[i].l;

				double ss = 0;
				for(int y = 0; y < 3; y++)
				{
					ss += ((2*TMath::Pi()/s.l)*pow(epsim, 3/2)*vol/3)*(dataJC_cor[i].e2/pow(L[y],2))/( pow(dataJC_cor[i].e1 + epsim*(1-L[y])/L[y], 2) + pow(dataJC_cor[i].e2, 2));
				}

				s.sig = ss;
				extJC.push_back(s);
				//printf("extJC = %f,%f \n",s.l, s.sig);
			}

			//calcolo del fattore di riscalamento per JC (best)
			double meanJC = 0;
			int nJC = 0;
			for (int i = 0; i < extJC.size(); ++i)
			{
				if(extJC[i].l >= 400 && extJC[i].l <= 420)
				{
				  meanJC += extJC[i].sig;
				  nJC++;
				}
			}
			meanJC /= nJC;

			double factor = meanSpet/meanJC;
			//cout << factor << endl;

			/*
			std::string label;
			ostringstream convert;   
			convert << r << "-" << epsim;
			label = convert.str();
			TString t_label = label;

			
			ofstream outf(t_label+"Avar.dat");
			for(int i =0 ; i < extJC.size() ; i++)
			{
				outf << extJC[i].l << " " << extJC[i].sig*factor<< endl;
			}
			outf.close();
			*/
			int m =0;
			for(int i =0 ; i < extJC.size() ; i++)
			{
				if(extJC[i].l >= 420 && extJC[i].l <= 650)
				{
					double expect = extJC[i].sig*factor;
					entry.chi2 += pow(expect - assSpet[i].sig,2)/expect;
					m++;
				}
			}
			chiTable.push_back(entry);
		 }
	}

	ofstream outf("chi2.csv");
	for(int i =0 ; i < chiTable.size() ; i++)
	{
		outf << chiTable[i].a1 << "," << chiTable[i].ecc << "," << chiTable[i].chi2 << endl;
	}
}
//   			cout << assSpet[i].l << " " << assSpet[i].sig << "---" << extJC[i].l << ":" << extJC[i].sig << << endl;
