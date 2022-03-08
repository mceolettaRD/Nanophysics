{

	struct diel {
		double l, E, e1, e2;
	};

   vector<diel> dataJC;

   int n=0;
   int m=0;
   
	ifstream inp;
	string fil="epsiJC.txt";
   inp.open(fil);

	double l,E,e1,e2;
	while(inp >> l >> E >> e1 >> e2)
	{
		diel d;
		d.l = l;
		d.E = E;
		d.e1 = e1;
		d.e2 = e2;
		dataJC.push_back(d);
		printf("%f:%f: %f,%f\n",d.l,d.E,d.e1,d.e2);
	}
   TGraph *g1 = new TGraph();
   TGraph *g2 = new TGraph();

   for (int i = 0; i < dataJC.size(); ++i)
   {
	  if(dataJC[i].l < 467)
	  {
		 g1->SetPoint(n++,dataJC[i].l, dataJC[i].e1);
		 g2->SetPoint(m++,dataJC[i].E, dataJC[i].e1);
	  }
   }
   g1->Fit("pol1","M");
   g1->Draw();

}