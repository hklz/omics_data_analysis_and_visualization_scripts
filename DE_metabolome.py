import os
import pandas as pd
import openpyxl
import regex as re_lib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#读取路径下所有.xlsx文件名字
def list_path_files(directory_path, suffix='.xlsx'):
    # Check if the directory exists
    if not os.path.isdir(directory_path):
        raise ValueError("路径下没这个文件")
    # List all files in the directory
    files = os.listdir(directory_path)
    # Filter out files that are not .xlsx
    files_find = [file for file in files if file.endswith(suffix)]
    return files_find

#转换xlsx文件到dataform
def read_excel_to_dataframe(file_path, columns=["geneID", "Regulation","GeneName",]):
    # Check if the file exists
    if not os.path.isfile(file_path):
        raise ValueError("The specified file does not exist.")
    # Read the Excel file
    df = pd.read_excel(file_path, usecols=columns, engine='openpyxl')
    return df

def count_elements_in_class(df, col='Class I'):
    element_counts = df[col].value_counts()
    return element_counts

def combine_counts(dataframes):
    """
    Combines the counts of elements in "Class I" across multiple DataFrames.

    Parameters:
    dataframes (list of pd.DataFrame): A list of DataFrames to process.

    Returns:
    pd.DataFrame: A DataFrame containing the combined counts of each element across all provided DataFrames.
    """
    combined_counts = pd.DataFrame()

    for i, df in enumerate(dataframes):
        # Count elements in the current DataFrame
        counts = df['Class I'].value_counts()

        # Convert to DataFrame and transpose
        counts_df = counts.to_frame().transpose()

        # Set the name of the index (which will become the column name after transposing)
        counts_df.index = [f'DataFrame_{i+1}']

        # Combine with the overall DataFrame
        combined_counts = pd.concat([combined_counts, counts_df], axis=0, sort=False)

    # Fill NaN values with 0, as they indicate no occurrence of that element in a DataFrame
    combined_counts = combined_counts.fillna(0).astype(int)

    return combined_counts

# 读取列表
DE_metabolome = "E:\sigma因子的研究\组学数据\衍生处理数据\DE_metabolome\\"
Output_path = "E:\sigma因子的研究\组学数据\衍生处理数据\DE_metabolome\\"
xlsx_files = list_path_files(DE_metabolome,".xlsx")
#print("输出读取的名字：")
#print(xlsx_files)

dfs = []

for file_name in xlsx_files:
    #轮流读取
    dataframe = read_excel_to_dataframe(DE_metabolome+file_name, columns=['Class I','Class II'])
    dfs.append(dataframe)
    print(file_name)
    #print(count_elements_in_class(dataframe,'Class I'))

com = combine_counts(dfs)
print(com)
com.to_excel(Output_path+'metabolome.xlsx', index=True)



