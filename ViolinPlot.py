import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 读取数据
df = pd.read_excel('E:/sigma因子的研究/组学数据/衍生处理数据/转录组热图聚类TPM.xlsx')

# 仅选择数值类型的列
df_numeric = df.select_dtypes(include=[np.number])

# 将数据转换为对数尺度，对数变换仅应用于数值列
df_log = np.log10(df_numeric + 1)  # 加1是为了避免对数0的情况,log10
#df_log = np.log(df_numeric + 1)  # 加1是为了避免对数0的情况,ln

# 重塑数据为长格式
df_long = df_log.reset_index().melt(id_vars='index', var_name='Sample', value_name='Value')


# 计算每个样本的描述性统计数据
stats = df_long.groupby('Sample')['Value'].describe()
print("Descriptive Statistics:\n", stats)
# 将统计数据保存到文件
stats.to_csv("descriptive_statistics.csv")

# 指定每个箱线图的颜色
palette = {
    "MGI": "#CDDC39",
    "MGD": "#FF9800",
    "MGE": "#2196F3",
    "MGF": "#F44336",
    "MGH": "#9C27B0",
    "MGN": "#4CAF50",
    "MGS": "#E91E63",
    "MGT": "#607D8B"
}
flierprops = dict(marker='o', color='black', markersize=2)

# 绘制箱线图
#sns.stripplot(data=df_log, size=2, color="black", jitter=True, alpha=0.5)
#sns.boxplot(x='Sample', y='Value', data=df_long, palette=palette,flierprops=flierprops)
sns.violinplot(x='Sample', y='Value', data=df_long, hue='Sample', palette=palette, legend=False)

plt.xticks(rotation=0)  # 旋转X轴标签以便清晰显示
plt.ylabel(r'$\log_{10}(\mathrm{TPM})$')
#plt.ylabel(r'$\ln(\mathrm{FPKM})$')
#plt.show()
plt.savefig(f"violinplot_TPM_log10", dpi=300, transparent=1)
plt.close()

