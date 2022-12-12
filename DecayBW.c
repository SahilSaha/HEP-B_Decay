Double_t myFunc ( Double_t x )
{
	Double_t gd= 0.002;
	Double_t E = 2.00697;
	Double_t pi = 3.1416;
	Double_t gamma ,k , denom ;
	gamma = sqrt ( x * x *( x * x + gd * gd ));
	k = 2* sqrt (2)* x * gd * gamma /( pi * sqrt ( x * x + gamma ));
	denom = 1/(( E *E - x * x )*( E *E - x * x )+ x * x * gd * gd );
	return k * denom ;
};


void DecayBW()
{
	gSystem->Load("libPhysics");

	TCanvas *c1 = new TCanvas();
	TCanvas *c2 = new TCanvas();
	TCanvas *c3 = new TCanvas();
	
	TLorentzVector B(0.0,0.0,0,5.27934);

	TF1 * BW = new TF1 ("BW"," myFunc (x)" ,0 ,2);


	TH1F *histed = new TH1F("histed","Energy-D", 100, 2, 3.3);
	TH1F *histep = new TH1F("histep","Energy-Pi", 100, 2, 3);
	TH1F *histpd = new TH1F("histpd","Momentum-D", 100, 2, 3);
	TH1F *histpp = new TH1F("histpp","Momentum-Pi", 100, 2, 3);
	
	TH1F *dec11 = new TH1F("dec11", "Momentum-D0", 100,2,2.2);
	TH1F *dec12 = new TH1F("dec12", "Momentum-Pi0", 100,0,0.3);
	TH1F *dec13 = new TH1F("dec13", "Energy-D0", 100,2.7,2.9);
	TH1F *dec14 = new TH1F("dec14", "Energy-Pi0", 100,0,0.5);
	
	TH1F *dec21 = new TH1F("dec21", "Momentum-D0", 100,1.6,2.5);
	TH1F *dec22 = new TH1F("dec22", "Momentum-gamma", 100,-0.2,0.5);
	TH1F *dec23 = new TH1F("dec23", "Energy-D0", 100,2.6,3.2);
	TH1F *dec24 = new TH1F("dec24", "Energy-gamma", 100,-0.2,0.5);
	
	for (Int_t n = 0; n < 1e6; n++)
	{
		Double_t m = BW->GetRandom();
		Double_t masses[2] = {m, 0.139570};

		TGenPhaseSpace event;
		event.SetDecay(B, 2, masses);
		Double_t weight = event.Generate();

		TLorentzVector *D = event.GetDecay(0);
		TLorentzVector *Pi = event.GetDecay(1);

		histed->Fill(D->Energy());
		histep->Fill(Pi->Energy());
		histpd->Fill(D->P());
		histpp->Fill(Pi->P());
		
		Double_t masses2[2] = {1.86484, 0.1349768};
		Double_t masses3[2] = {1.86484, 0.0};

		TGenPhaseSpace event1;
		TGenPhaseSpace event2;
		event1.SetDecay(*D, 2, masses2);	
		
		Double_t weight1 = event1.Generate();

		TLorentzVector *D01 = event1.GetDecay(0);
		TLorentzVector *Pi0 = event1.GetDecay(1);	
		
		event2.SetDecay(*D, 2, masses3);	
		
		Double_t weight2 = event2.Generate();

		TLorentzVector *D02 = event2.GetDecay(0);
		TLorentzVector *gamma = event2.GetDecay(1);
		
		dec11->Fill(D01->P());
		dec12->Fill(Pi0->P());
		dec13->Fill(D01->E());
		dec14->Fill(Pi0->E());
		
		dec21->Fill(D02->P());
		dec22->Fill(gamma->P());
		dec23->Fill(D02->E());
		dec24->Fill(gamma->E());
	}
	
	c1->Divide(2,2);
	c1->cd(1); histed->Draw();
	c1->cd(2); histep->Draw();
	c1->cd(3); histpd->Draw();
	c1->cd(4); histpp->Draw();
	
	c2->Divide(2,2);
	c2->cd(1); dec11->Draw();
	c2->cd(2); dec12->Draw();
	c2->cd(3); dec13->Draw();
	c2->cd(4); dec14->Draw();
	
	c3->Divide(2,2);
	c3->cd(1); dec21->Draw();
	c3->cd(2); dec22->Draw();
	c3->cd(3); dec23->Draw();
	c3->cd(4); dec24->Draw();
}	
	

