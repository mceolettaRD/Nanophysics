{
	const double z = 1; //lunghezza cuvetta
	struct assp {
	double l, E, ass;
	};

	vector<assp> absorb;

	ifstream inp;
	string fil="AuNP.txt";
   inp.open(fil);

	double l,ass;
   	while(inp >> l >> ass)
   	{
   		assp d;
   		d.l = l;
   		d.E = 1239.84193/l;
   		d.ass = ass;
   		absorb.push_back(d);
   		printf("%f:%f,%f\n",d.l,d.E,d.ass);
   	}

   	ofstream outf("sigmaExt_prod_AuNP.txt");

   	for(int i = 0; i < absorb.size(); i++)
   	{
   		double prod = absorb[i].ass/(z*log(TMath::E())/log(10)) ;
   		printf("%f:%f, %f\n",absorb[i].l, absorb[i].E, prod);
   		outf << absorb[i].l << " " << absorb[i].E << " " << prod << endl;
   	}

   	outf.close();
}