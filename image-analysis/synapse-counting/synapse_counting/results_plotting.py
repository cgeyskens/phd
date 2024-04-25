import pandas as pd
import seaborn as sns

def data_formatting(df, id_vars, value_vars, value_name, var_name = "condition"):
        """
        Wrapper function for formating the data for plotting into seaborn swarmplot.

        Args:
           df: the dataframe from which to plot
           id_vars: column with the image filename containing text for plotting on the x-axis (eg CA1 SR, CA1 SLM, ...)
           value_vars: columns to unpivot (see documentation pd.melt)
           value_name: the name on the y-axis
           var_name: the condition of rotated vs actual in synaptic marker colocalization

        Returns:
           df_melted: a dataframe pivoted and ready as input into seaborn swarmplot.
        """
        df_melted = pd.melt(df, id_vars = id_vars, value_vars = value_vars, var_name = var_name, value_name = value_name)
        df_melted["condition"] = df_melted[var_name].apply(lambda x: "rotated" if value_vars[1] in x else "actual")
        df_melted["hippocampal layer"] = df_melted[id_vars[0]].apply(lambda x: " ".join(x.split("_")[-2:]))
        return df_melted

def plot_data(df, x, y, extra_y, title):  
    p = sns.swarmplot(x=x, y=y, hue = "condition",
                    data = df,
                    # jitter = False,
                    dodge = True,
                    marker = "o",
                    alpha = 0.5)
    max_y = df[y].max()
    y_upper_limit = max_y + extra_y
    p.set_title(title)
    p.set_ylim(0, y_upper_limit)
    return p
