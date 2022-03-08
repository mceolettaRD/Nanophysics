{
	struct assp {
	double l, E, ass;
	};

   const double lmin = 515;
   const double lmax = 540;

   double n = 0;

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

      TGraph *spettro = new TGraph();

      for (int i = 0; i < absorb.size(); ++i)
      {
         if(absorb[i].l > lmin && absorb[i].l < lmax)
         {
            spettro->SetPoint(n++, absorb[i].l, absorb[i].ass);
         }
      }

      spettro->SetMarkerStyle(3);
      spettro->Draw("AP");
      TF1 *fun =  new TF1("fun","gaus");
      spettro->Fit(fun,"M");
}