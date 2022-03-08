{
	const double r = 10;
	const double epsim = 1.77;
	const double vol = 4*TMath::Pi()*pow(r,3)/3;

	struct diel {
		double l, E, e1, e2;
	};

	struct sigma {
	double l, sig;
	};

	int n=0;

   	vector<diel> dataJC;

	ifstream inp;
	string fil="epsiJC.txt";
   	inp.open(fil);

	double l,e1,e2,E;
   	while(inp >> l >> E >> e1 >> e2)
   	{
   		diel d;
   		d.l = l;
   		d.E = 1239.84193/l;
   		d.e1 = e1;
   		d.e2 = e2;
   		dataJC.push_back(d);
   		printf("%f:%f,%f\n",d.l,d.e1,d.e2);
   	}

   	vector<sigma> ext;

   	ofstream outf("sigmaExtJC.txt");

   	for(int i = 0; i < dataJC.size(); i++)
   	{
   		sigma s;
   		s.l = dataJC[i].l;
   		s.sig = 9*(2*TMath::Pi()/s.l)*pow(epsim, 3/2)*vol*(dataJC[i].e2/(pow(dataJC[i].e1 + 2*epsim, 2) + pow(dataJC[i].e2, 2)));
   		ext.push_back(s);
   		printf("%f:%f, %f\n",s.l, dataJC[i].E, s.sig);
   		outf << s.l << " " << dataJC[i].E << " " << s.sig << endl;
   	}

   	outf.close();
}