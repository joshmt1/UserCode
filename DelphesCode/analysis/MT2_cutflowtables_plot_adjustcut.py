# goal
# for each lumi, choose the best cut based on Zbi
# make a plot that shows the Zbi as a function of lumi with a 'dynamic' cut

#args: <NMn where n=1,2,3> <PhaseI/PhaseII>

#event counts for 300fb-1 are hard-coded below

# == instructions for making summary plot ==
##run this code 3 times:
#python MT2_cutflowtables_plot_adjustcut.py NM1
#python MT2_cutflowtables_plot_adjustcut.py NM2
#python MT2_cutflowtables_plot_adjustcut.py NM3
## this will general 3 root files with TGraphs in a subdirectory called DelphesMT2 (must exist before running script)
## then run MT2_cutflowtables_plot_adjustcut_savepdf_combined.C to generate final output

import ROOT,sys
from ROOT import TGraph,TLegend,TCanvas,RooStats,TFile
cutval = [400,500,600,700,800]
#phaseII numbers
n_sm_2  = [3.53,1.59,0.75,0.38,0.23]
n_nm1_2 = [14.1,9.9,6.4,4.3,2.8]
n_nm2_2 = [46.1,27.7,16.2,9.1,5.0]
n_nm3_2 = [53.1,40.8,31.7,21.9,13.6]

#phaseI numbers
n_sm_1  = [ 4.61 , 2.23,  1.36,  1.01, 0.39 ]
n_nm1_1 = [ 15.3 , 10.8,  7.1,  4.4,  3.0 ]
n_nm2_1 = [ 49.1 , 29.4 , 17.2 , 9.9 , 5.5 ]
n_nm3_1 = [ 59.0 , 46.3 , 35.1 , 23.5 , 15.2 ]

if len(sys.argv)<3:
    print 'give me arguments!'
    sys.exit(1)
which = sys.argv[1]
if not 'NM' in which:
    sys.exit(1)
phase=sys.argv[2]

if phase=='PhaseI':
    print 'PhaseI'
    n_sm=n_sm_1
    n_nm1=n_nm1_1
    n_nm2=n_nm2_1
    n_nm3=n_nm3_1
else:
    print 'PhaseII'
    n_sm=n_sm_2
    n_nm1=n_nm1_2
    n_nm2=n_nm2_2
    n_nm3=n_nm3_2

#update: we want to use only the 500,600,800 points
#these were seen to mostly provide the best sensitivity while minimizing the number of search regions
#cutval = [ cutval[1], cutval[2],cutval[4]]
#n_sm = [ n_sm[1], n_sm[2],n_sm[4]]
#n_nm1 = [ n_nm1[1], n_nm1[2],n_nm1[4]]
#n_nm2 = [ n_nm2[1], n_nm2[2],n_nm2[4]]
#n_nm3 = [ n_nm3[1], n_nm3[2],n_nm3[4]]


def zbi(b,s,f):
    return RooStats.NumberCountingUtils.BinomialObsZ(b+s,b,f)

def get_graph( sm300, sig300, db):
    g1 = TGraph()
    g2 = TGraph()
    
    for ip,lumi in enumerate(range(10,3001,10)):
        scale = float(lumi)/float(300)
        bestfom = -99
        bestcut = -99
        for ism,isig,icut in zip(sm300,sig300,cutval):
            fom = zbi(scale*ism,scale*isig,db)
            if fom>bestfom:
                bestfom=fom
                bestcut=icut
#        print lumi,bestcut,bestfom
        g1.SetPoint(ip,lumi,bestfom)
        g2.SetPoint(ip,lumi,bestcut)
    return (g1,g2)



if which=='NM1':
    S = n_nm1
if which=='NM2':
    S = n_nm2
if which=='NM3':
    S = n_nm3


(g_550,cut_550) = get_graph(n_sm, S, 0.3)
g_550.SetName('mt2_550')
cut_550.SetName('mt2_550_bestcut')

#sys.exit(1)

(g_550_0p15,cut_550_0p15) = get_graph(n_sm, S, 0.15)
g_550_0p15.SetName('mt2_550_0p15')
cut_550_0p15.SetName('mt2_550_0p15_bestcut')

(g_550_0p50,cut_550_0p50) = get_graph(n_sm, S, 0.5)
g_550_0p50.SetName('mt2_550_0p50')
cut_550_0p50.SetName('mt2_550_0p50_bestcut')

## plot comparing cut values
#g_550.SetLineColor(ROOT.kRed)
#g_550.Draw('AL')

#g_500.SetLineColor(ROOT.kGreen)
#g_500.Draw('L')

#g_450.SetLineColor(ROOT.kBlue)
#g_450.Draw('L')


## plot comparing systematics values
thecanvas=TCanvas('thecanvas','',600,600)

#do not draw 15%
#g_550_0p15.SetLineColor(ROOT.kBlue)
#g_550_0p15.Draw('al')

linewidth = 3

g_550.SetLineColor(ROOT.kBlue)
g_550.SetLineWidth(linewidth)
g_550.Draw('al')

g_550_0p50.SetLineColor(ROOT.kRed)
g_550_0p50.SetLineWidth(linewidth)
g_550_0p50.Draw('l')

g_550.GetHistogram().SetXTitle("Luminosity (fb^{-1})")
g_550.GetHistogram().SetYTitle("Zbi")
g_550.GetHistogram().SetMinimum(0)

leg_x1 = 0.47
leg_x2 = 0.87
leg_y1 = 0.24
leg_y2 = 0.5

leg=TLegend(leg_x1,leg_y1,leg_x2,leg_y2)
#leg.AddEntry(g_550_0p15,'15% Systematics','l')
legtext = '#splitline{%d%s systematic}{uncertainty}' 
leg.AddEntry(g_550,legtext % (30,'%'),'l')
leg.AddEntry(g_550_0p50,legtext % (50,'%'),'l')
leg.SetBorderSize(0);
leg.SetLineStyle(0);
leg.SetTextFont(42);
leg.SetFillStyle(0);
leg.Draw()

#filename='DelphesMT2/MT2_significance_'+which
#filename+='.pdf'
#thecanvas.SaveAs(filename)

#plot best cut value
thecanvas2=TCanvas('thecanvas2','',600,600)

cut_550_0p50.Draw('al')
cut_550_0p50.SetLineWidth(linewidth)
cut_550_0p50.SetLineColor(ROOT.kRed)
cut_550_0p50.GetHistogram().SetMinimum(300)

cut_550.SetLineColor(ROOT.kBlue)
cut_550.SetLineWidth(linewidth)
cut_550.Draw('l')

#cut_550_0p15.SetLineColor(ROOT.kBlue)
#cut_550_0p15.Draw('l')

cut_550_0p50.GetHistogram().SetXTitle("Luminosity (fb^{-1})")
cut_550_0p50.GetHistogram().SetYTitle("MT2 cut value")

leg2=TLegend(leg_x1,leg_y1,leg_x2,leg_y2)
#leg.AddEntry(cut_550_0p15,'15% Systematics','l')
leg2.AddEntry(cut_550,legtext % (30,'%'),'l')
leg2.AddEntry(cut_550_0p50,legtext % (50,'%'),'l')
leg2.SetBorderSize(0);
leg2.SetLineStyle(0);
leg2.SetTextFont(42);
leg2.SetFillStyle(0);
leg2.Draw()

#filename='DelphesMT2/MT2_significance_bestcut_'+which
#filename+='.pdf'
#thecanvas2.SaveAs(filename)

which+=phase
filename='DelphesMT2/MT2_significance_'+which+'.root'
outfile = TFile(filename,'recreate')
g_550_0p50.Write()
thecanvas.Write()
thecanvas2.Write()
outfile.Close()

