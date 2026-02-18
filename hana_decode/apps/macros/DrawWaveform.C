void DrawWaveform()
{
  gErrorIgnoreLevel = kWarning;
  Int_t runnumber = 242;

  const Int_t mode = 10;
  const UInt_t slot = 3;
  const UInt_t NUMSAMPLE = 240 / 4;
  const UInt_t ns0 = 40;
  const UInt_t ns1 = 120;
  const UInt_t ng = 50;
  const UInt_t threshold = 600;

  TGraph *g_sample[8];
  for (Int_t ig = 0; ig < 8; ig++)
    g_sample[ig] = new TGraph(ng);

  TH1 *h_spectrum[8];
  for (Int_t ih = 0; ih < 8; ih++)
    h_spectrum[ih] = new TH1F(Form("h_spectrum_%d", ih), Form("ADC spectrum for PMT %d", ih + 1), 90, 10.5, 100.5);

  const char *type[3] = {"Left", "Right", "All"};
  TH1 *h_sum[3];
  for (Int_t ih = 0; ih < 3; ih++)
    h_sum[ih] = new TH1F(Form("h_sum_%d", ih), Form("ADC spectrum for %s PMT Sum", type[ih]), 290, 10.5, 300.5);

  auto f = new TFile(Form("Rootfiles/fadc_data_%d.root", runnumber));
  TDirectory *dir = (TDirectory *)f->Get(Form("/mode_%d_data/slot_%u", mode, slot));
  TTree *t_store = (TTree *)dir->Get("waveform");
  UInt_t store_event, store_channel, store_sample[100];
  t_store->SetBranchAddress("event", &store_event);
  t_store->SetBranchAddress("channel", &store_channel);
  t_store->SetBranchAddress("sample", store_sample);

  UInt_t last_event = 0;
  UInt_t total_channel = 0;
  UInt_t fadc_channel = 0;
  vector<UInt_t> v_channel;
  UInt_t max_index[8] = {};
  UInt_t max_sample[8] = {};
  UInt_t max_sample_1 = 0;
  UInt_t max_sample_2 = 0;
  UInt_t event_1or2 = 0;
  UInt_t event_1and2 = 0;
  Float_t sum_sample[8] = {};
  bool trig[2] = {};

  TString wavefile = Form("Waveform-run%d.pdf", runnumber);
  auto c0 = new TCanvas("c0", "c0", 600, 600);
  c0->Print(wavefile + "[");

  for (ULong64_t ien = 0; ien < t_store->GetEntries(); ien++)
  {
    t_store->GetEntry(ien);

    if (store_event != last_event)
    {
      last_event = store_event;
      ien--;

      if (trig[0] && trig[1])
      {
        for (ULong64_t jen = ien + 1 - total_channel; jen <= ien; jen++)
        {
          t_store->GetEntry(jen);
          // cout << store_event << ", " << store_channel << ", " << store_sample[0] << endl;
          if (store_channel < 8)
          {
            v_channel.push_back(store_channel);
            for (UInt_t sample_num = 0; sample_num < ng; sample_num++)
              g_sample[store_channel]->SetPoint((Int_t)sample_num, sample_num * 4, store_sample[sample_num]);
          } // PMT channels
        } // jen

        c0->cd();
        auto leg0 = new TLegend(0.1, 0.65, 0.25, 0.9);
        for (UInt_t ich = 0; ich < v_channel.size(); ich++)
        {
          UInt_t chan = v_channel.at(ich);
          g_sample[chan]->SetTitle(Form("Event %u", store_event));
          g_sample[chan]->GetXaxis()->SetTitle("Time (ns)");
          g_sample[chan]->GetYaxis()->SetRangeUser(0, 1600);
          g_sample[chan]->SetLineColor(chan + 1);
          g_sample[chan]->SetLineStyle(1);
          g_sample[chan]->SetLineWidth(3);
          g_sample[chan]->Draw(ich == 0 ? "AL" : "L");
          leg0->AddEntry(g_sample[chan], Form("CH%u", chan), "L");
          leg0->Draw();
        }
        c0->Print(wavefile);
        c0->Clear("D");
        delete leg0;
        v_channel.clear();

        for (Int_t ic = 0; ic < 8; ic++)
          h_spectrum[ic]->Fill(sum_sample[ic]);

        UInt_t sum_left = 0;
        for (Int_t ic = 0; ic < 4; ic++)
          sum_left += sum_sample[ic];
        h_sum[0]->Fill(sum_left);

        UInt_t sum_right = 0;
        for (Int_t ic = 4; ic < 8; ic++)
          sum_right += sum_sample[ic];
        h_sum[1]->Fill(sum_right);
        h_sum[2]->Fill(sum_left + sum_right);
      } // trig[0] && trig[1]

      if (max_sample_1 > threshold || max_sample_2 > threshold)
        event_1or2++;
      if (max_sample_1 > threshold && max_sample_2 > threshold)
        event_1and2++;

      total_channel = -1;
      fadc_channel = 0;
      for (Int_t ic = 0; ic < 8; ic++)
      {
        max_index[ic] = 0;
        max_sample[ic] = 0;
        max_sample_1 = 0;
        max_sample_2 = 0;
        sum_sample[ic] = 0;
      }
      for (Int_t it = 0; it < 2; it++)
        trig[it] = 0;
    } // new event

    else if (store_channel < 8)
    {
      for (Int_t is = 0; is < NUMSAMPLE; is++)
      {
        if (store_sample[is] > max_sample[store_channel])
        {
          max_index[store_channel] = is;
          max_sample[store_channel] = store_sample[is];
        }
        if (is >= ns0 / 4 && is < ns1 / 4)
          sum_sample[store_channel] += store_sample[is];
      }
      for (Int_t is = 40 / 4; is < 120 / 4; is++)
        if (store_sample[is] > max_sample_1)
          max_sample_1 = store_sample[is];
      for (Int_t is = 120 / 4; is < 200 / 4; is++)
        if (store_sample[is] > max_sample_2)
          max_sample_2 = store_sample[is];

      const Int_t nped = 4;
      Float_t ped = 0.;
      for (Int_t is = 0; is < nped; is++)
        ped += store_sample[is];
      max_sample[store_channel] -= (Float_t)ped / nped;
      sum_sample[store_channel] /= (Float_t)(ns1 - ns0) / 4;
      sum_sample[store_channel] -= (Float_t)ped / nped;
      sum_sample[store_channel] *= (Float_t)(ns1 - ns0) / 4 / NUMSAMPLE;

      if (max_sample[store_channel] > threshold)
      {
        fadc_channel++;
        trig[store_channel / 4] = true;
      }
    } // PMT channels

    total_channel++;
  } // ien

  c0->Print(wavefile + "]");

  cout << "OR event = " << event_1or2 << "; AND event = " << event_1and2 << endl;

  TString specfile = Form("Spectrum-run%d-%dto%dns.pdf", runnumber, ns0, ns1);
  auto c1 = new TCanvas("c1", "c1", 4 * 600, 2 * 600);
  c1->Divide(4, 2);
  for (Int_t ih = 0; ih < 8; ih++)
  {
    c1->cd(ih + 1);
    h_spectrum[ih]->Draw("HIST");
  }
  c1->Print(specfile + "(");

  auto c2 = new TCanvas("c2", "c2", 2 * 600, 2 * 600);
  c2->Divide(2, 2);
  for (Int_t ih = 0; ih < 3; ih++)
  {
    c2->cd(ih + 1);
    gPad->SetLogy();
    h_sum[ih]->Draw("HIST");
  }
  c2->Print(specfile + ")");
}
