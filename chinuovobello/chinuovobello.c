{

	ifstream inp;		//leggo i dati
   	double x,y,sx,sy;
   	double xv[10000];
   	double yv[10000];

   	double res[10000];
   	int n=0;

   	string fil="input.dat";
   	cout << "File:\n";

   	inp.open(fil);

   	cout.precision(15);

   	while(inp >> x >> y)
   	{
   		xv[n]=x;
   		yv[n]=y;
   		cout << n << "\n"<<endl;
		cout << "Letti  x: "<< x <<" y:"<< y << "\n";
   		n++;
   	}

   	inp.close();


	h = new TGraph(n, xv, yv);		
	h->SetTitle("fit; X; Counts");

	for(int i=0; i<n;i++) 
	{
		h->SetPoint(i, xv[i], yv[i]);
	}
	
	f = new TF1("f", "2*[0]/([1]*sqrt((4*[2]/([3]*[1]*[1]))-1))*exp(-x*[1]/(2*[2]))*sin(x*[1]*(1/(2*[2]))*sqrt((4*[2]/([3]*[1]*[1]))-1)+[4])");	//funzione ha n parametri, ma solo [1] [2] [3] vengono steppati sono gli unici
	f->SetParameters(21,60,1e-3,100e-9,-2); //guess dei parametri (il fit deve convergere...)

	h->Fit(f,"M"); 

	h->SetMarkerStyle(3);

	R = f->GetParameter(1);
	L = f->GetParameter(2);
	C = f->GetParameter(3);

	cout<<"parametri best fit steppati [1] [2] [3] (RLC):  " << R <<" : "<< L << " : " << C <<endl;
	
	double chi2n=0;
	for(int i=0; i<n; i++)
	{
		double A=f->Eval(xv[i]);
		double res2n=pow((yv[i]-A),2);
		chi2n+=res2n;
		//cout<<"predetta"<< f->Eval(xv[i])<< "  osservata :" << yv[i]<< "  +- " <<syv[i]<<endl;

	}

	int gdl=f->GetNDF();
	cout <<"gradi di libertà : "<<gdl<<endl;

	double chi=f->GetChisquare() ;
	cout<<"  chi2 : "<<chi2n<<"   rootchi2:"<<chi<<endl;
 
	

double step1,step2,stepg;
cout<<"step traslazione [3] [2] [1] (CLR, [1] è graficoso):  "<<endl; //le variazioni: step 1 è di C, step 2 è di L, stepg è di R
cin >> step1 >> step2>> stepg;

int m;
cout<<"numero iterazioni normali simmetriche:"<<endl; //punti generati nello spazio dei parametri non graficosi
cin>>m;

int u;
cout<<"numero iterazioni graficose simmetriche (se vuoi tenere [1] al best fit fanne 0):"<<endl; //file creati, best fit centrale e gli altri attorno (2u+1)
cin>>u;

double e_a1=-1*(step1*m); //valori di partenza per griglia sia graficosa che normale
double e_a2=-1*(step2*m);
double e_ag=-1*(stepg*u); 




double R0=R-stepg*u; //parto dal minimo valore graficoso e lo imposto, poi verrà variato
f->SetParameter(1,R0);


for(int k=-1*u; k<=u;k++) //ciclo graficoso, ogni iterazione produce un file di output per un dato valore di R fisso
{

	vector <double> zchi; //vettori del chi2 e degli altri due parametri
	vector <double> xnu; 
	vector <double> ynu;

	for (int i=-1*m;i<=m;i++) //ciclo su C
	{
		e_a1+=step1; //incremento il valore dell'offset
		f->SetParameter(3,C+e_a1); //imposto il parametro ofssettato C nella funzione di fit
		f->GetParameter(3); //controllo che sia giusto
		for (int j=-1*m;j<=m;j++)
		{
			e_a2+=step2; //idem per  L
			f->SetParameter(2,L+e_a2);
			f->GetParameter(2);
			xnu.push_back(f->GetParameter(3));		//assex C mi scrivo nei vettori i parametri
			ynu.push_back(f->GetParameter(2));		//assey L

			double chi2m=0; //calcolo del chi2 e lo inserisco nel vettore
			
			for(int i=0; i<n; i++)
			{	
				double A=f->Eval(xv[i]);
				double res2m=pow((yv[i]-A),2);
				chi2m+=res2m;
			}
			zchi.push_back(chi2m);
		}
		f->SetParameter(2,L); //reimposto L al valore di prima perchè ora ho finito un ciclo su L e quindi devo cambiare C di nuovo e farne un'altro sempre dall'inizio
		e_a2=-1*(step2*m); //anche l'offset di L deve tornare all'inizio
	}

stringstream ss;
ss << R0;
string str = ss.str();

string nomefile = str+"chi"+fil;

ofstream outputfile(nomefile); //scrivo i vettori in un file

for (int i=0;i<pow(2*m+1,2);i++){
outputfile<<xnu[i]<<", " << ynu[i] << ", " << zchi[i]<<endl;
}

outputfile.close();

f->SetParameter(3,C); //rimetto la funzione nelle ci con i parametri di best fit
f->SetParameter(2,L);
R0+=stepg;	//faccio lo stepping sul parametro graficoso e lo reimposto
e_a1=-1*(step1*m);
f->SetParameter(1,R0);
cout<< "iterazione graficosa k = " << k*stepg << endl;
} //fine ciclo graficoso

}








































