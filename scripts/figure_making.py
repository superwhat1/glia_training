# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 15:58:46 2024

@author: BioCraze
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scypi.stats as stats


#Figures for non-viral transfection project

gcamp_summary = 
transfected_count = 

fig,axes =plt.subplots(nrows=2, ncols=3, figsize = (25,10), layout="constrained")
palette = ["#95beff","#00bac6"]
sns.set_palette(palette)

ax1 = sns.barplot(x = "transfection_technique", y = "fluo per tectal slice", hue = "Trial", data=gcamp_summary, capsize=0.2, ax=axes[0,0])
ax1.set_ylim(top=max(gcamp_summary["fluo per tectal slice"])*1.2)
sns.stripplot(x = "transfection_technique", y = "fluo per tectal slice", data = gcamp_summary, color="k", size=5, hue="Trial", dodge=True, legend=False, ax=ax1)
ax2=sns.barplot(x = "transfection_technique", y = "fluorescence per ROI", hue = "Trial", data=gcamp_summary, capsize=0.2, ax=axes[0,1])
ax2.set_ylim(top=max(gcamp_summary["fluorescence per ROI"])*1.2)
sns.stripplot(x = "transfection_technique", y = "fluorescence per ROI", data = gcamp_summary, color="k", size=5, hue="Trial", dodge=True, legend=False, ax=ax2)
ax3=sns.barplot(x = "transfection_technique", y = "cell to slice fluo ratio", hue = "Trial", data=gcamp_summary, capsize=0.2, ax=axes[0,2])
ax3.set_ylim(top=max(gcamp_summary["cell to slice fluo ratio"])*1.2)
sns.stripplot(x = "transfection_technique", y = "cell to slice fluo ratio", data = gcamp_summary, color="k", size=5, hue="Trial", dodge=True, legend=False, ax=ax3)
ax4=sns.barplot(x = "transfection_technique", y = "cell count", hue = "Trial", data=gcamp_summary, capsize=0.2, ax=axes[1,0])
ax4.set_ylim(top=max(gcamp_summary["cell count"])*1.2)
sns.stripplot(x = "transfection_technique", y = "cell count", data = gcamp_summary, color="k", size=5, hue="Trial", dodge=True, legend=False, ax=ax4)
ax5=sns.barplot(x = "transfection_technique", y = "per cell fluo", hue = "Trial", data=cell_fluo, capsize=0.2, ax=axes[1,1])
#ax5.set_ylim(top=max(cell_fluo["per cell fluo"])*1.2)
#sns.stripplot(x = "transfection_technique", y = "per cell fluo", data = cell_fluo, color="k", size=3, hue="Trial", dodge =True,ax=ax5)
ax6=sns.barplot(x = "Transfection technique", y = "transfected cell count", data=transfected_count, capsize=0.2, palette=["#00a258","#827000","#ad0f00"], ax=axes[1,2])
ax6.set_ylim(top=max(transfected_count["transfected cell count"])*1.2)
sns.stripplot(x = "Transfection technique", y = "transfected cell count", data = transfected_count, color="k", size=5, legend=False, ax=ax6)
plt.show()
#------------------------------------------------------------------------------------------

#Figures for glia plasticity
#GLIA DATA
    #no training
fig,axes =plt.subplots(nrows=2, ncols=4, figsize = (25,10), layout="constrained")
palette = ["#95beff","#00bac6"]
sns.set_palette(palette)

responsive_glia_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_glia_responding_cell_numbers.csv").sort_values(by=['time group']).replace(to_replace='min99',value='min100')
responsive_glia_data=responsive_glia_data.loc[~((responsive_glia_data['time group']=='min30') | (pd.isna(responsive_glia_data['time group'])))]
responsive_glia_ctrl = sns.pointplot(x = "time group", y = "Corrected cell count", hue = "treatment", data=responsive_glia_data[responsive_glia_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[0,0])
responsive_glia_ctrl.set_ylim(top=max(responsive_glia_data["Corrected cell count"])*1.2, bottom=min(responsive_glia_data["Corrected cell count"])*0.6)
responsive_glia_ctrl.set_title("Active glia count - no training")
sns.stripplot(x = "time group", y = "Corrected cell count", data = responsive_glia_data[responsive_glia_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=responsive_glia_ctrl)

responsive_glia_exp= sns.pointplot(x = "time group", y = "Corrected cell count", hue = "treatment", data=responsive_glia_data[responsive_glia_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[1,0])
responsive_glia_exp.set_ylim(top=max(responsive_glia_data["Corrected cell count"])*1.2, bottom=min(responsive_glia_data["Corrected cell count"])*0.6)
responsive_glia_exp.set_title("Active glia count - training")
sns.stripplot(x = "time group", y = "Corrected cell count", data = responsive_glia_data[responsive_glia_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=responsive_glia_exp)


glia_transients_count_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_glia_response_numbers.csv").sort_values(by=['time group']).replace(to_replace='min99',value='min100')
glia_transients_count_data=glia_transients_count_data.loc[~((glia_transients_count_data['time group']=='min30') | (pd.isna(glia_transients_count_data['time group'])))]
glia_transients_count_ctrl = sns.pointplot(x = "time group", y = "normalized response count", hue = "treatment", data=glia_transients_count_data[glia_transients_count_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[0,1])
glia_transients_count_ctrl.set_ylim(top=max(glia_transients_count_data["normalized response count"])*1.2, bottom=min(glia_transients_count_data["normalized response count"])*0.6)
glia_transients_count_ctrl.set_title("Glia transients counts - no training")
sns.stripplot(x = "time group", y = "normalized response count", data = glia_transients_count_data[glia_transients_count_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_count_ctrl)
	
	#With training
glia_transients_count_exp= sns.pointplot(x = "time group", y = "normalized response count", hue = "treatment", data=glia_transients_count_data[glia_transients_count_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[1,1])
glia_transients_count_exp.set_ylim(top=max(glia_transients_count_data["normalized response count"])*1.2, bottom=min(glia_transients_count_data["normalized response count"])*0.6)
glia_transients_count_exp.set_title("Glia transients counts - training")
sns.stripplot(x = "time group", y = "normalized response count", data = glia_transients_count_data[glia_transients_count_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_count_exp)

glia_transients_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_glia_summary_by_cell.csv").sort_values(by=['time group']).replace(to_replace='min99',value='min100')
glia_transients_data=glia_transients_data.loc[~((glia_transients_data['time group']=='min30') | (pd.isna(glia_transients_data['time group'])))]
#glia_transients_amp_ctrl = sns.pointplot(x = "time group", y = "peaks", hue = "treatment", data=glia_transients_data[glia_transients_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[0,2])
#glia_transients_amp_ctrl.set_ylim(top=max(glia_transients_data["peaks"])*1.2, bottom=min(glia_transients_data["peaks"])*0.6)
#glia_transients_amp_ctrl.set_title("Glia peak amplitude - no training")
#sns.stripplot(x = "time group", y = "peaks", data = glia_transients_data[glia_transients_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_amp_ctrl)

glia_transients_amp_exp= sns.pointplot(x = "time group", y = "peaks", hue = "treatment", data=glia_transients_data[glia_transients_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[1,2])
glia_transients_amp_exp.set_ylim(top=max(glia_transients_data["peaks"])*1.2, bottom=min(glia_transients_data["peaks"])-min(glia_transients_data["peaks"])*1.2)
glia_transients_amp_exp.set_title("Glia peak amplitude - training")
sns.stripplot(x = "time group", y = "peaks", data = glia_transients_data[glia_transients_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_amp_exp)

glia_transients_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_glia_summary_by_cell.csv").sort_values(by=['time group']).replace(to_replace='min99',value='min100')
glia_transients_data=glia_transients_data.loc[~((glia_transients_data['time group']=='min30') | (pd.isna(glia_transients_data['time group'])))]
#glia_transients_amp_ctrl = sns.pointplot(x = "time group", y = "area", hue = "treatment", data=glia_transients_data[glia_transients_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[0,3])
#glia_transients_amp_ctrl.set_ylim(top=max(glia_transients_data["area"])*1.2, bottom=min(glia_transients_data["area"])*0.6)
#glia_transients_amp_ctrl.set_title("Glia auc - no training")
#sns.stripplot(x = "time group", y = "area", data = glia_transients_data[glia_transients_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_amp_ctrl)

glia_transients_amp_exp= sns.pointplot(x = "time group", y = "area", hue = "treatment", data=glia_transients_data[glia_transients_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[1,3])
glia_transients_amp_exp.set_ylim(top=max(glia_transients_data["area"])*1.2, bottom=min(glia_transients_data["area"])-min(glia_transients_data["area"])*1.2)
glia_transients_amp_exp.set_title("Glia auc - training")
sns.stripplot(x = "time group", y = "area", data = glia_transients_data[glia_transients_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=glia_transients_amp_exp)

plt.show()


#NEURON DATA
    #No training
fig,axes =plt.subplots(nrows=2, ncols=1, figsize = (10,10), layout="constrained")
palette = ["#e51bec","#1be6ec"]
sns.set_palette(palette)

neuron_response_data= pd.read_csv("E:/glia projects/plasticity/summaries/glia_training_neuron_summary_by_cell.csv").sort_values(by=['time group']).replace(to_replace='min99',value='min100')
neuron_response_data=neuron_response_data.loc[~((neuron_response_data['time group']=='min30') | (pd.isna(neuron_response_data['time group'])))]
neuron_response_amp_ctrl = sns.pointplot(x = "time group", y = "peaks", hue = "treatment", data=neuron_response_data[neuron_response_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[0])
neuron_response_amp_ctrl.set_ylim(top=max(neuron_response_data["peaks"])*1.2, bottom=min(neuron_response_data["peaks"])*0.6)
neuron_response_amp_ctrl.set_title("Neuron peak amplitude - no training")
sns.stripplot(x = "time group", y = "peaks", data = neuron_response_data[neuron_response_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=neuron_response_amp_ctrl)

neuron_response_auc_ctrl = sns.pointplot(x = "time group", y = "area", hue = "treatment", data=neuron_response_data[neuron_response_data['treatment'].str.contains('without training')], errorbar='se', capsize=0.2, ax=axes[1])
neuron_response_auc_ctrl.set_ylim(top=max(neuron_response_data["area"])*1.2, bottom=min(neuron_response_data["area"])*0.6)
neuron_response_auc_ctrl.set_title("Neuron auc - no training")
sns.stripplot(x = "time group", y = "area", data = neuron_response_data[neuron_response_data['treatment'].str.contains('without training')], size=5, hue="treatment", dodge=False, legend=False, ax=neuron_response_auc_ctrl)
plt.show()
    
    #With training
fig,axes =plt.subplots(nrows=2, ncols=1, figsize = (10,10), layout="constrained")
palette = ["#e51bec","#1be6ec"]
sns.set_palette(palette)

neuron_response_amp_exp= sns.pointplot(x = "time group", y = "peaks", hue = "treatment", data=neuron_response_data[neuron_response_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[0])
neuron_response_amp_exp.set_ylim(top=max(neuron_response_data["peaks"])*1.2, bottom=min(neuron_response_data["peaks"])-min(neuron_response_data["peaks"])*1.2)
neuron_response_amp_exp.set_title("Neuron peak amplitude - training")
sns.stripplot(x = "time group", y = "peaks", data = neuron_response_data[neuron_response_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=neuron_response_amp_exp)

neuron_response_auc_exp= sns.pointplot(x = "time group", y = "area", hue = "treatment", data=neuron_response_data[neuron_response_data['treatment'].str.contains('with training')], errorbar='se', capsize=0.2, ax=axes[1])
neuron_response_auc_exp.set_ylim(top=max(neuron_response_data["area"])*1.2, bottom=min(neuron_response_data["area"])-min(neuron_response_data["area"])*1.2)
neuron_response_auc_exp.set_title("Neuron auc - training")
sns.stripplot(x = "time group", y = "area", data = neuron_response_data[neuron_response_data['treatment'].str.contains('with training')], size=5, hue="treatment", dodge=False, legend=False, ax=neuron_response_auc_exp)
plt.show()



#CELL LOCATION MATCHING
fig,ax =plt.subplots(figsize = (25,25), layout="constrained")
ax = sns.heatmap(fltrd_dst_matrix, vmin=np.min(fltrd_dst_matrix), vmax=np.min(fltrd_dst_matrix)*100 ,center=np.min(fltrd_dst_matrix)*10, cbar=False)
ax.set_ylabel("cell id in current time point")
ax.set_xlabel("cell id in next time point")
ax.set_title("Matching cells")

#FIELD POTENTIAL GRAPHS
NBQX_APV = pd.read_csv("E:/glia projects/field recordings/field_potential_data.csv")

fig,ax =plt.subplots(figsize = (25,10), layout="constrained")
palette = ["#95beff","#00bac6"]
sns.set_palette(palette)

ax = sns.pointplot(x = "series time", y = "Normalized amp", hue = "treatment", data=NBQX_APV, errorbar='se', capsize=0.2, ax=ax)
ax.set_ylim(top=max(NBQX_APV["Normalized amp"])*1.2, bottom=min(NBQX_APV["Normalized amp"])-min(NBQX_APV["Normalized amp"])*1.2)
ax.set_title("Field recordings")
sns.stripplot(x = "series time", y = "Normalized amp", data = NBQX_APV, size=5, hue="treatment", dodge=False, legend=False, ax=ax)

plt.show()



