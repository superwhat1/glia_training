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

