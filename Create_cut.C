void run10555cuts361(){
  TFile *fout=new TFile("dis_cut_361.root","recreate");//create rootfile on which you will save your cuts
  TCutG *cutg1 = new TCutG("mycut_ph_dp",11);
  cutg1->SetVarX("dp");
  cutg1->SetVarY("ph");
  cutg1->SetPoint(0,-0.0388423, 0.0326031);
  cutg1->SetPoint(1,-0.045172, 0.0168999);
  cutg1->SetPoint(2,-0.0477912, -0.0129868);  
  cutg1->SetPoint(3,-0.0432076, -0.0347686);
  cutg1->SetPoint(4,-0.0233456, -0.0355284);
  cutg1->SetPoint(5,0.0526104, -0.0246375);
  cutg1->SetPoint(6,0.0567575, -0.0124802);
  cutg1->SetPoint(7,0.0561027, 0.0113278);
  cutg1->SetPoint(8,0.0486817, 0.0196859);
  cutg1->SetPoint(9,0.0179064, 0.0265244);
  cutg1->SetPoint(10,-0.0388423, 0.0326031);
  
  TCanvas * c1 = new TCanvas("c1", "c1");
  cutg1->Draw();  
  cutg1->Write();//save first cut
  
  
  
  
  TCutG *cutg2 = new TCutG("mycut_ph_th",16);
  cutg2->SetVarX("th");
  cutg2->SetVarY("ph");    
  cutg2->SetPoint(0,-0.0519382, 0.0162206);
  cutg2->SetPoint(1,-0.0508469, 0.00732186);
  cutg2->SetPoint(2,-0.0602322, 0.00470457);
  cutg2->SetPoint(3,-0.0602322, -0.00707322);
  cutg2->SetPoint(4,-0.0549939, -0.0253942);
  cutg2->SetPoint(5,-0.0471364, -0.0327226);
  cutg2->SetPoint(6,0.0526104, -0.0319374);
  cutg2->SetPoint(7,0.0619958, -0.024609);
  cutg2->SetPoint(8,0.0678889, -0.00838186);
  cutg2->SetPoint(9,0.0665794, 0.0049663);
  cutg2->SetPoint(10,0.0606862, 0.00391939);
  cutg2->SetPoint(11,0.0591584, 0.0264281);
  cutg2->SetPoint(12,0.0523922, 0.0285219);
  cutg2->SetPoint(13,-0.0425528,0.0266898);
  cutg2->SetPoint(14,-0.0499738, 0.0222404);
  cutg2->SetPoint(15,-0.0519382, 0.0162206);
  
  TCanvas * c2 = new TCanvas("c2", "c2");
  cutg2->Draw();
  cutg2->Write();//save second cut




  TCutG *cutg3 = new TCutG("mycut_ph_y",11);
  cutg3->SetVarX("y");
  cutg3->SetVarY("ph");    
  cutg3->SetPoint(0,-0.0246551, 0.0331096);
  cutg3->SetPoint(1,-0.0272743, 0.027031);
  cutg3->SetPoint(2,-0.0272743, -0.0208383);
  cutg3->SetPoint(3,-0.0240003, -0.0251441);
  cutg3->SetPoint(4,0.0248909, -0.0367948);
  cutg3->SetPoint(5,0.0310023, -0.0355284);
  cutg3->SetPoint(6,0.0347128, -0.0248908);
  cutg3->SetPoint(7,0.0377685, 0.014873);
  cutg3->SetPoint(8,0.0364589, 0.0196859);
  cutg3->SetPoint(9,-0.0152698, 0.0333629);
  cutg3->SetPoint(10,-0.0246551, 0.0331096);
 
  TCanvas * c3 = new TCanvas("c3", "c3");
  cutg3->Draw();
  cutg3->Write();//save 3rd cut


  TCutG *cutg4 = new TCutG("mycut_th_dp",10);
  cutg4->SetVarX("dp");
  cutg4->SetVarY("th");    
  cutg4->SetPoint(0,-0.0473546, 0.0650225);
  cutg4->SetPoint(1,-0.0381875, -0.0471791);
  cutg4->SetPoint(2,-0.0248734, -0.0573102);
  cutg4->SetPoint(3,0.0488999, -0.0568036);
  cutg4->SetPoint(4,0.0565392, -0.0509783);
  cutg4->SetPoint(5,0.0563209, -0.0332489);
  cutg4->SetPoint(6,0.0316571, 0.0541316);
  cutg4->SetPoint(7,0.00721145, 0.0665421); 
  cutg4->SetPoint(8,-0.0423346, 0.0670487); 
  cutg4->SetPoint(9,-0.0473546, 0.0650225);  
  TCanvas * c4 = new TCanvas("c4", "c4");
  cutg4->Draw(); 
  cutg4->Write();//save 4th cut



  TCutG *cutg5 = new TCutG("mycut_th_y",9);
  cutg5->SetVarX("y");//th
  cutg5->SetVarY("th"); //y   
  cutg5->SetPoint(0,-0.0274926, 0.0604528);
  cutg5->SetPoint(1,-0.0259647, -0.0560165);
  cutg5->SetPoint(2,0.034058, -0.0562782);
  cutg5->SetPoint(3,0.0360224, -0.0376955);
  cutg5->SetPoint(4,0.0358041, 0.0473664);
  cutg5->SetPoint(5,0.0329666, 0.064117);
  cutg5->SetPoint(6,0.0299109, 0.0664726);
  cutg5->SetPoint(7,-0.0218177, 0.0672577);
  cutg5->SetPoint(8,-0.0274926, 0.0604528);
  
  TCanvas * c5 = new TCanvas("c5", "c5");
  cutg5->Draw(); 
  cutg5->Write();//save 5th cut



  TCutG *cutg6 = new TCutG("mycut_y_dp",9);
  cutg6->SetVarX("dp");//dp
  cutg6->SetVarY("y"); //y  
  cutg6->SetPoint(0,-0.0445172, 0.0356424);
  cutg6->SetPoint(1,-0.0466998, 0.000943456);
  cutg6->SetPoint(2,-0.0436441, -0.0266637);
  cutg6->SetPoint(3,0.0523922, -0.0253973);
  cutg6->SetPoint(4,0.0558844, -0.00589502);
  cutg6->SetPoint(5,0.0534835, 0.0333629);
  cutg6->SetPoint(6,0.0416972, 0.0376686);
  cutg6->SetPoint(7,-0.0403702, 0.0369088);
  cutg6->SetPoint(8,-0.0445172, 0.0356424);
  
  TCanvas * c6 = new TCanvas("c6", "c6");
  cutg6->Draw();
  cutg6->Write();//save 6th cut
  fout->Close();//now close your output rootfile
  
}
