import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import numpy as np
import io
from PIL import Image
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas


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


def plot_data(df, x, y, extra_y_upper, title, hue = "gRNA", extra_y_lower = 0, ax = None):
    """
    Wrapper function for custom swarmplot from seaborn.

    Args:
        df: the dataframe outputted from the data_formatting function.
        x: data for x-axis
        y: the metric analyzed
        extra_y_upper: how much extra space above the highest y value in all conditions
        title: title of the plot
        hue: in this case always gRNA
        extra_y_lower: how much extra space below the lowest y value in all conditions
    
    Returns:
        p: a custom swarmplot.
    """
    if ax is None:
        fig, ax = plt.subplots()
    p = sns.stripplot(x = x, 
                      y = y, 
                      hue = hue,
                      data = df,
                      dodge = True,
                      marker = "o",
                      alpha = 0.5,
                      ax = ax)
    max_y = df[y].max()
    min_y = df[y].min()

    y_upper_limit = max_y + extra_y_upper
    
    if min_y < 0:
        y_lower_limit = min_y + extra_y_lower
    else:
        y_lower_limit = 0
    
    p.set_title(title)
    p.set_ylim(y_lower_limit, y_upper_limit)

    return p


class PlotResults:

    def __init__(self, df, candidate_gRNA, name_of_plot, control_gRNA = "LacZ-gRNA"):
        self.df = df
        self.candidate_gRNA = candidate_gRNA    
        self.control_gRNA = control_gRNA
        self.name_of_plot = name_of_plot

        """
        Args:
            df: the dataframe with all the data, including the brain nr
            candidate_gRNA: which is the candidate_gRNA from which you want to check the statistics
            control_gRNA: is standard LacZ-gRNA
            name_of_plot: name of the plot provided
        """

    def check_statistics(self, hippocampal_layer, metric):
        """
        This function checks the different between hemispheres per brain. Uses a paired t-test to compare the means.

        Args:
            df: the dataframe with all the data
            hippocampal_layer: which hippocampal layer to check the statistics from
            metric: whichc metric to check the statistics from
            candidate_gRNA: which is the candidate_gRNA from which you want to check the statistics
            control_gRNA: is standard LacZ-gRNA
        
        Returns:
            filtered_df_hip_layer: the filtered dataframe with the hippocampal layer
            replicate_averages_long: dataframe with the means calculated from each brain, in compact format
            replicate_averages_wide: dataframe with the means calculated from each brain, in long format
            statistic: the T statistic of the paired t-test
            pvalue: the accompanying p-value of the paired t-test
        """
        filtered_df_hip_layer = self.df[self.df["hippocampal_layer"].str.contains(hippocampal_layer)]
        replicate_averages_long = filtered_df_hip_layer.groupby(["gRNA", "Brain"], as_index = False).agg({metric:"mean"})
        replicate_averages_wide = replicate_averages_long.pivot_table(columns = "gRNA", values = metric, index = "Brain")
        statistic, pvalue = scipy.stats.ttest_rel(replicate_averages_wide[self.control_gRNA], replicate_averages_wide[self.candidate_gRNA])
        
        return filtered_df_hip_layer, replicate_averages_long, replicate_averages_wide, statistic, pvalue
    

    def create_scatter_plot_with_means(self, ax, df_averages, df_points, metric):
        """"
        Creates a single scatter plot (one metric & one hippocampal layer), where the brain means are plotted and connected with lines.
        Points on the side represent each individual measured data point.
        
        Args:
            ax: no need to fullfill argument, is used for plotting the mean and individual points.
            df_averages: dataframe with the means calculated from each brain, in long format (relates to the "replicate_averages_long" output from check_statistics).
            df_points: dataframe that contains the individual points per brain from a specific hippocampal layer (related to the "filtered_df_hip_layer" output from check statistics).
            metric: which metric to plot.
            candidate_gRNA: which is the candidate_gRNA from which you want to check the statistics.
            control_gRNA: is standard LacZ-gRNA.

        Returns:
            A single scatter plot where one metric is plotted for one hippocampal layer. 
        """
        # individual points - preparing the data for plotting
        df_points_LacZ_gRNA = df_points[df_points["gRNA"] == self.control_gRNA]
        df_points_cand_gRNA = df_points[df_points["gRNA"] == self.candidate_gRNA]

        # mean points - preparing the data for plotting
        df_average_LacZ_gRNA = df_averages[df_averages["gRNA"] == self.control_gRNA]
        df_average_cand_gRNA = df_averages[df_averages["gRNA"] == self.candidate_gRNA]

        # set the individual points across the x axis
        x_LacZ_points = np.ones(len(df_points_LacZ_gRNA)) * 0.1  
        x_VCAM1_points = np.ones(len(df_points_cand_gRNA)) * 0.9  

        # set the mean points across the x axis
        x_LacZ_averages = np.ones(len(df_average_LacZ_gRNA)) * 0.25  
        x_VCAM1_averages = np.ones(len(df_average_cand_gRNA)) * 0.75 

        # plot the individual points
        ax.scatter(x_LacZ_points, df_points_LacZ_gRNA[metric], color = "#808080", edgecolors = "black", alpha = 0.6, linewidths = 0.5)
        ax.scatter(x_VCAM1_points, df_points_cand_gRNA[metric], color = "#c92ffb", edgecolors = "black", alpha = 0.6, linewidths = 0.5)

        # plot the mean points
        ax.scatter(x_LacZ_averages, df_average_LacZ_gRNA[metric], color = "#808080", edgecolors = "black", linewidths = 0.8)
        ax.scatter(x_VCAM1_averages, df_average_cand_gRNA[metric], color = "#c92ffb", edgecolors = "black", linewidths = 0.8)

        # plot the lines connecting the mean points
        for i in range(len(df_average_LacZ_gRNA)):
            ax.plot([0.25, 0.75], [df_average_LacZ_gRNA.iloc[i][metric], df_average_cand_gRNA.iloc[i][metric]], linewidth = 0.5, c = "k")

        # other attributes
        ax.spines[['right', 'top']].set_visible(False)


    def create_scatter_plot_with_means_per_hippocampal_layer(self, hippocampal_layer_list, metric, ax = None):
        """
        Creates a scatter plot (one metric & all hippocampal layers), where the brain means are plotted and connected with lines. 
        Points on the side represent each individual measured data point.
        Also plots the p value above the each comparison and a star if significant. 
        
        Args:
            df: the large dataframe containing all the metric measurement from all the hippocampal layers, needs to have a column that specifies which brains is measured.
            hippocampal_layer_list: list of hippocampal layers that needs to be plotted.
            metric: which metric to plot.
            candidate_gRNA: which is the candidate_gRNA from which you want to check the statistics.
            control_gRNA: is standard LacZ-gRNA.
        
        Returns:
            Scatter plots where one metric is plotted for all hippocampal layers. 
        """

        # set the size of the figure
        if ax is None:
            fig, axes = plt.subplots(nrows = 1, ncols = len(hippocampal_layer_list), figsize=(12, 6))

        # variables to track the global y-axis limits
        global_y_max = float('-inf')

        for idx, hippocampal_layer in enumerate(hippocampal_layer_list):
            
            # get statistics
            df_filtered, df_averages, _, _, p_value = self.check_statistics(hippocampal_layer = hippocampal_layer, metric = metric)


            # update global y-axis limits
            local_y_max = max(df_filtered[metric].max(), df_averages[metric].max())
            global_y_max = max(global_y_max, local_y_max)
            
            # make the plot on the specific subplot axis
            ax

            # make the plot on the specific subplot axis
            ax = axes[idx]
            self.create_scatter_plot_with_means(ax, df_averages = df_averages, df_points = df_filtered, metric = metric)

            # set x-axis labels to hippocampal layer names
            ax.set_xlabel(f"{hippocampal_layer}")

            # remove x-tick labels
            ax.set_xticks([])

            # remove y-axis labels and y axis lines for all but the first plot
            if idx > 0:
                ax.set_ylabel("")
                ax.set_yticks([])
                ax.spines['left'].set_visible(False)

            # add p-value to each graph
            p_value_str = f"p = {p_value:.3f}"
            p_value_text = ax.text(0.5, 0.95, p_value_str, ha = 'center', va='bottom', transform=ax.transAxes, fontsize=12)  # Font size for p-value

            # add significant stars as a separate text object with a larger font size
            stars = ""
            if p_value < 0.0001:
                stars = "****"
            elif p_value < 0.001:
                stars = "***"
            elif p_value < 0.01:
                stars = "**"
            elif p_value < 0.05:
                stars = "*"
            stars_text = ax.text(0.5, 0.98, stars, ha='center', va='bottom', transform=ax.transAxes, fontsize=20)  # Larger font size for stars

        # apply global y-axis limits to all subplots
        global_y_max = global_y_max * 1.1  # adjust the maximum y-axis to be 0.1 higher
        for ax in axes:
            ax.set_ylim(0, global_y_max)

        # create a color-coded legend outside the loop with patch objects for color matching
        handles = [
            plt.Rectangle((0, 0), 1, 1, color="#808080", ec="black", lw=0.5),
            plt.Rectangle((0, 0), 1, 1, color="#c92ffb", ec="black", lw=0.5),
        ]
        labels = [self.control_gRNA, self.candidate_gRNA]
        fig.legend(handles, labels, loc='upper left', title="gRNA", bbox_to_anchor=(1, 1))  # Adjust x and y coordinates as needed

        # Add a title to the figure
        fig.suptitle(metric, fontsize=18)

        plt.tight_layout()

        return fig
    

    def save_plot_to_bytes(self, fig):
        """
        Helper function that saves an image into buffer memory. 
        This is way around the axes from matplotlib because its already used to create the figure from the
        "create_scatter_plot_with_means_per_hippocampal_layer" function.

        Args:
            fig: scatter plots where one metric is plotted for all hippocampal layers (relates to the output from "create_scatter_plot_with_means_per_hippocampal_layer" )
        
        Returns:
            image of the object.
        """
        buf = io.BytesIO()
        canvas = FigureCanvas(fig)
        fig.savefig(buf, format='png', bbox_inches='tight') # to ensure that the legend for each plot is printed out.
        buf.seek(0)
        return Image.open(buf)
    

    def save_figure(self, hippocampal_layer_list, metric_list):
        """
        This is the master function, where it all comes together to save the plot.

        Args:
            hippocampal_layer_list: list of hippocampal layers assessed
            metric_list: list of metrics assessed.
        """
        # loops over the images
        images = []
        for metric in metric_list:
            fig = self.create_scatter_plot_with_means_per_hippocampal_layer(hippocampal_layer_list = hippocampal_layer_list, metric = metric)
            image = self.save_plot_to_bytes(fig)
            images.append(image)
            plt.close(fig) # do not print out each figure

        # extracts figure dimensions
        widths, heights = zip(*(img.size for img in images))
        total_height = sum(heights)
        max_width = max(widths) +1

        # combines the figures and pastes it into the blank canvas
        combined_image = Image.new('RGB', (max_width, total_height))
        y_offset = 0
        for img in images:
            combined_image.paste(img, (0, y_offset))
            y_offset += img.height

        # saves the image
        combined_image.save(self.name_of_plot + ".png")