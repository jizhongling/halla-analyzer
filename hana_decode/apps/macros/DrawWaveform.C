void DrawWaveform(const Int_t runnumber)
{
  gErrorIgnoreLevel = kWarning;

  const Int_t mode = 1;
  const UInt_t slot = 3;
  const UInt_t ng = 100;

  TGraph *g_sample = new TGraph(ng);

  auto f = new TFile(Form("Rootfiles/fadc_data_%d.root", runnumber));
  TDirectory *dir = (TDirectory *)f->Get(Form("/mode_%d_data/slot_%u", mode, slot));
  TTree *t_store = (TTree *)dir->Get("waveform");
  UInt_t store_event, store_channel, store_sample[1000];
  t_store->SetBranchAddress("event", &store_event);
  t_store->SetBranchAddress("channel", &store_channel);
  t_store->SetBranchAddress("sample", store_sample);

  TString wavefile = Form("plots/Waveform-run%d.pdf", runnumber);
  auto c0 = new TCanvas("c0", "c0", 600, 600);
  c0->Print(wavefile + "[");

  for (ULong64_t ien = 0; ien < t_store->GetEntries(); ien++)
  {
    t_store->GetEntry(ien);

    //cout << "Event " << store_event << ", Chan " << store_channel << ": ";
    for (UInt_t sample_num = 0; sample_num < ng; sample_num++)
    {
      //cout << store_sample[sample_num] << " ";
      g_sample->SetPoint((Int_t)sample_num, sample_num * 4, store_sample[sample_num]);
    }
    //cout << endl;

    c0->cd();
    g_sample->SetTitle(Form("Event %u, Chan %u", store_event, store_channel));
    g_sample->GetXaxis()->SetTitle("Time (ns)");
    g_sample->SetLineStyle(1);
    g_sample->SetLineWidth(3);
    g_sample->Draw("AL");
    c0->Print(wavefile);
    c0->Clear("D");
  } // ien

  c0->Print(wavefile + "]");
}
