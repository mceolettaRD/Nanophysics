{
	//const double epsim = 2.182; //fisso valore di epsim a quello di frolich
	const double e0 = 9.84;
	const double gammab = 0.0654926;
	const double Ep = 9.010;
	const double vf = 1.40e15; //in nm/s
	const double hbar = 6.582119569e-16; //eV/s

	const double rho = 2;
	const double r = 8.27;

	struct diel {
		double l, E, e1, e2;
	};

	struct sigma {
	double l, sig;
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


			double vol = 4*TMath::Pi()*pow(r,3)/3; //volume NP
			//printf("%f, %f \n",r,epsim);

			double gammar = gammab + vf/(4*r)*hbar;

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
				s.sig = 9*(2*TMath::Pi()/s.l)*pow(epsim, 3/2)*vol*(dataJC_cor[i].e2/(pow(dataJC_cor[i].e1 + 2*epsim, 2) + pow(dataJC_cor[i].e2, 2)));
				extJC.push_back(s);
				//printf("extJC = %f,%f \n",s.l, s.sig);
			}

			//calcolo del fattore di riscalamento per JC (best)
			double factor = 1e7*rho*log10(TMath::E());
			//cout << factor << endl;

			
			std::string label;
			ostringstream convert;   
			convert << r << "-" << rho;
			label = convert.str();
			TString t_label = label;

			
			ofstream outf(t_label+"Uniq.dat");
			for(int i =0 ; i < extJC.size() ; i++)
			{
				outf << extJC[i].l << " " << extJC[i].sig*factor<< endl;
			}
			outf.close();
			

}
//   			cout << assSpet[i].l << " " << assSpet[i].sig << "---" << extJC[i].l << ":" << extJC[i].sig << << endl;
