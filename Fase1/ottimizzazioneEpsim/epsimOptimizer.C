{

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


   for (double epsim = 1.67; epsim <= 2.5; epsim += 0.001)
   {
      for(int i = 0; i < dataJC.size(); i++)
      {
         if(dataJC[i].l == 527)
         {
            cout << epsim << ":" << dataJC[i].e1 + 2*epsim << endl;
         }
      }
   }



}