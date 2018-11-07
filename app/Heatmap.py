from app.settings import *
import math, os, json, csv, argparse
import pandas as pd
from collections import defaultdict, Counter
from textwrap import wrap
import matplotlib # this needs to be added to run on galaxylab
matplotlib.use('Agg') # this needs to be added to run on galaxylab
from matplotlib import gridspec
import seaborn as sns
import matplotlib.pyplot as plt

class Heatmap(object):
    """
    This is a program that genreates a heatmap of multiple RGI analyses.
    """

    def __init__(self, input, classification, frequency, output, cluster, display, debug):
        self.input = input
        self.classification = classification
        self.frequency = frequency
        self.output = output
        self.cluster = cluster
        self.display = display
        self.debug = debug

        if self.debug:
            logger.setLevel(10)

    def __repr__(self):
        """Returns Heatmap class full object."""
        return "Heatmap({}".format(self.__dict__)

    def get_figure_dimensions(self,s,g):
        """Set the dimensions of the figure"""
        w = len(s)
        l = len(g)
        figsize = (w, l)
        fig = plt.figure(figsize = figsize)
        # print(figsize)
        return w,l,fig,figsize

    def create_plot(self,t,r): # t = type plot, r = ratio
        # ax0 = heatmap, ax1 = categories, ax2 = frequency
        """Creates the appropriate number of subplots"""
        if t == 'c': # Category
            gs = gridspec.GridSpec(1, 2, width_ratios=[1,r])
            ax0 = plt.subplot(gs[1]) # Heatmap
            ax1 = plt.subplot(gs[0], sharey=ax0) # Categories
            ax1.set_xlim([0,1])
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)
            ax1.tick_params(bottom=False, left=False)
            return ax0,ax1,gs
        if t == 'cf': # Category and frequency
            gs = gridspec.GridSpec(2, 2, width_ratios=[1,r], height_ratios=[1,50])
            ax0 = plt.subplot(gs[3])
            ax1 = plt.subplot(gs[2], sharey=ax0)
            ax1.set_xlim([0,1])
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.spines['bottom'].set_visible(False)
            ax1.spines['left'].set_visible(False)
            ax1.tick_params(bottom=False, left=False)
            ax2 = plt.subplot(gs[1], sharex=ax0)
            ax2.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.tick_params(bottom=False)
            plt.setp(ax2.get_xticklabels(), visible=False)
            ax2.set_axisbelow(True)
            plt.grid(axis='y')
            return ax0,ax1,ax2,gs
        if t == 'f': # Frequency
            gs = gridspec.GridSpec(2, 1, height_ratios=[1,50])
            ax0 = plt.subplot(gs[1])
            ax2 = plt.subplot(gs[0], sharex=ax0)
            ax2.set_xlim([0,1])
            ax2.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.tick_params(bottom=False)
            plt.setp(ax2.get_xticklabels(), visible=False)
            return ax0,ax2,gs

    def create_class_series(self, class_dict, name):
        """Create a pandas series for the classification chosen"""
        class_df = pd.Series(class_dict, name=name)
        class_df.index.name = "model_name"
        class_df = class_df.apply(tuple)
        class_df.reset_index()
        return class_df

    def create_categories(self, class_dict, df):
        """Reformats the dataframe to handle categorization data"""
        for model in class_dict:
            if len(class_dict[model]) > 1:
                df = df.append([df.loc[model]]*(len(class_dict[model])-1))

        # Assigns a unique identifier to each entry to index the dataframe without duplicates
        count = Counter(df.index.values)
        new_labels = df.index.tolist() # Add * to the models with duplicates
        new_index = []
        counted = {}
        for model in list(df.index.values):
            if count[model] > 1:
                idx = new_labels.index(model)
                new_labels[idx] = "%s*" % (model)

        for i,v in enumerate(list(df.index.values)):
            if v in counted:
                counted[v] += 1
                new_index.append(v+"_"+str(counted[v]))
            else:
                counted[v] = 0
                new_index.append(v+"_0")

        df.index = new_labels
        df = df.assign(uID=new_index)
        df = df.reset_index().set_index("uID")
        return df

    def create_frequency_df(self, df, outfile):
        """Creates a dataframe for frequency data"""
        freq_df = pd.DataFrame() # New matrix df based on frequencies
        freq_dict = {} # Dictionary to keep track of ocurance of resistome profile
        samples = {} # Dictionary to group samples with identifcal profiles together
        n = 0
        for column in df:
            if column != 'index':
                n += 1
            s1 = df.loc[:, column] # Store column data as a Series
            if freq_df.empty:
                freq_df = pd.concat([freq_df, s1], axis = 1, sort=True)
                freq_dict[column] = 1
                samples[column] = [column]
            else:
                counter = 0
                for profile in freq_df:
                    # print(profile)
                    s2 = df.loc[:, profile]
                    if s1.equals(s2):
                        counter += 1
                        freq_dict[profile] += 1
                        samples[profile].append(column)
                        break
                if counter == 0:
                    freq_df = pd.concat([freq_df, s1], axis=1, sort=True)
                    freq_dict[column] = 1
                    samples[column] = [column]
        try:
            del freq_dict["index"]
            del samples["index"]
        except:
            pass

        # Order columns by frequency of resistome profiles
        cols = sorted(samples.keys(), key=(lambda x: len(samples[x])), reverse=True)
        if self.classification:
            cols.insert(0, 'index')
        freq_df = freq_df[cols]

        # Create .txt file of grouped samples
        with open('{o}-{n}-frequency.txt'.format(o=outfile, n=n), 'w') as f:
            fcsv = csv.writer(f, delimiter='\t')
            fcsv.writerow(['Frequency', 'Samples'])
            for s in cols:
                if s != 'index':
                    fcsv.writerow([len(samples[s])] + [', '.join(map(str, samples[s]))])
        return freq_df, freq_dict

    def draw_barplot(self, freq_dict, ax2, order):
        """Draws the frequency barplot"""
        from matplotlib.ticker import MaxNLocator
        y = list(freq_dict.values())
        yint = range(min(y), math.ceil(max(y)))
        order_values = [freq_dict[x] for x in order]
        bp = ax2.bar(range(len(freq_dict)), order_values, color="k", align="edge")
        ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax2.set_ylabel("Profile Frequency", rotation=0, va='center', labelpad=150, fontsize='xx-large')

    def cluster_data(self, option, df):
        """Hierarchically clusters the dataframe"""
        if option == "samples":
            cm = sns.clustermap(df, row_cluster=False, col_cluster=True)
            global clustered_col
            clustered_col = cm.dendrogram_col.reordered_ind
            df = df.iloc[:, clustered_col]
            # clustered_col is a list of indexes
        elif option == "genes":
            cm = sns.clustermap(df, row_cluster=True, col_cluster=False)
            clustered_row = cm.dendrogram_row.reordered_ind
            df = df.iloc[clustered_row, :]
        elif option =="both":
            cm = sns.clustermap(df, row_cluster=True, col_cluster=True)
            clustered_col = cm.dendrogram_col.reordered_ind
            clustered_row = cm.dendrogram_row.reordered_ind
            df = df.iloc[clustered_row, clustered_col]
        return df

    def calculate_categories(self, series, width):
        """Creates category labels"""
        freq = series.value_counts()
        freq = freq.sort_index()
        categories = freq.index.values
        # Introduces a line break if the category name is too long
        if float(width) < 7.0:
            categories = ['\n'.join(wrap(cat,40)) for cat in categories]
        elif 7.0 < float(width) < 7.5:
            categories = ['\n'.join(wrap(cat,50)) for cat in categories]
        else: # If category axes is larger, wrap less
            categories = ['\n'.join(wrap(cat,55)) for cat in categories]
        ranges = freq.values
        return categories, ranges

    def draw_categories(self, ax1, ranges, cat_list, ax0, display, df):
        """Draws all elements of categories on axes"""
        pal = sns.color_palette("dark") # for text
        light_pal = sns.color_palette("pastel") # for fill
        label_bb = ax0.yaxis.get_ticklabels()[0].get_window_extent() # get height of labels

        # Get longest gene name to determine width of marker fill
        longest_gene = max(df.index.tolist(), key=len) # string
        i = df.index.tolist().index(longest_gene) # index
        fill_width = ax0.yaxis.get_ticklabels()[i].get_window_extent().width

        # Draw first category first
        ax1.plot([1,1], [0,ranges[0]], lw=10,
            color=pal[3])
        ax1.text(0.5, (ranges[0]/2), cat_list[0], horizontalalignment="center",
            fontsize='xx-large', color=pal[3], weight="bold")

        if display == "fill":
            for line in ax0.yaxis.get_ticklines()[0:ranges[0]*2]:
                line.set_color(light_pal[3])
                line.set_markersize(fill_width)
                # print(label_bb.height)
                line.set_markeredgewidth(label_bb.height-(label_bb.height*0.25))
                # print(label_bb.height-(label_bb.height*0.05))
        elif display == "text":
            for label in ax0.yaxis.get_ticklabels()[0:ranges[0]]:
                label.set_color(pal[3])
        i = 0

        # Automate drawing of the rest of the categories
        c_iter = 0
        for x in ranges[1:]: #skips first item
            i += 1
            if c_iter > 3: # Reset iterator to cycle through 4 colours
                c_iter = 0
            ymax = int(math.fsum(ranges[0:i]) + x)
            ymin = int(math.fsum(ranges[0:i]))
            ax1.plot([1,1], [ymin, ymax], lw = 10, color=pal[c_iter])
            if cat_list[i].count('\n') > 1:
                temp_str = " ".join(cat_list[i].split('\n'))
                new_str = '\n'.join(wrap(temp_str,60))
                ax1.text(0.5, (ymin + (ymax - ymin)/2), new_str, fontsize='large', horizontalalignment="center", color=pal[c_iter], weight="bold")
            else:
                ax1.text(0.5, (ymin + (ymax - ymin)/2), cat_list[i], fontsize='xx-large', horizontalalignment="center", color=pal[c_iter], weight="bold")

            if display == "fill":
                for line in ax0.yaxis.get_ticklines()[ymin*2:ymax*2]:
                    line.set_clip_on(True)
                    line.set_color(light_pal[c_iter])
                    line.set_markersize(fill_width)
                    line.set_markeredgewidth(label_bb.height-(label_bb.height*0.25))
            elif display == "text":
                for label in ax0.yaxis.get_ticklabels()[ymin:ymax]:
                    label.set_color(pal[c_iter])

            c_iter += 1

    def get_axis_size(self, fig, ax):
        """Returns the width and length of a subplot axes"""
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        return width, height

    def write_csv(self, category, df, class_series, output):
        """Produces a tab-delimited file of heatmap matrix"""
        # Convert series to dataframe
        class_frame = class_series.to_frame()
        # Inner join between series and df on unique IDs, rename gene column
        merged_df = pd.merge(df, class_frame, left_index=True, right_index=True)
        merged_df.rename(columns={'index':'gene'}, inplace=True)
        # Moves classification column to the front
        col_list = merged_df.columns.tolist()
        if self.classification in col_list:
            reordered_columns = [self.classification] + [v for i,v in enumerate(col_list) if i != col_list.index(self.classification)]
        else:
            logger.warning("Couldn't create csv of heatmap matrix. " +
            "The category {} did not exist.".format(self.classification))
        merged_df = merged_df[reordered_columns]
        merged_df.to_csv('{}.csv'.format(output), index=False)

    def run(self):
        # print args
        logger.info(json.dumps(self.__dict__, indent=2))

        # List to hold the file name
        directory = self.input
        files = os.listdir(directory)
        jsons = []
        shortened = []
        for thing in files:
            file_path = os.path.join(directory, thing)
            if thing.endswith(".json") and os.path.isfile(file_path): # Check if it's a file
                jsons.append(thing)
        genelist = [] # List of unique genes
        genes = {} # Will become the dataframe
        resist_mech = {} # key: gene, value: resistance mechanism
        drug_class = {} # key: gene, value: drug class
        gene_family = {} # key: gene, value: gene family
        excluded = [] # incompletely curated models
        for jsonfile in jsons:
            # {json file: {Model: type_hit}}
            accession = jsonfile.split(".json")[0]
            shortened.append(accession) # Don't take whole file name
            genes[accession] = {}
            with open(os.path.join(directory, jsonfile)) as data: # Use os.path.join
                rgi_data = json.load(data)
            try:
                del rgi_data["_metadata"]
            except:
                pass

            try:
                tophits = {}
                # Top hit of each ORF
                for key,value in rgi_data.items():
                    if isinstance(value, dict):
                        contig_id = key
                        hsp = max(value.keys(), key=(lambda key: value[key]['bit_score']))

                        # Flag to exclude loose hits
                        if value[hsp]["type_match"] != "Loose":
                            topmodel = value[hsp]["model_name"]
                            tophits[topmodel] = value[hsp]["type_match"]

                            # Build dictionary of model names and their classifications
                            try:
                                if self.classification:
                                    rm = 0
                                    gf = 0
                                    dc = 0
                                    for entry in value[hsp]["ARO_category"]:
                                        if value[hsp]["ARO_category"][entry]["category_aro_class_name"] == "Resistance Mechanism":
                                            rm += 1
                                            if value[hsp]["model_name"] not in resist_mech:
                                                resist_mech[value[hsp]["model_name"]] = [value[hsp]["ARO_category"][entry]["category_aro_name"]]
                                            else:
                                                if value[hsp]["ARO_category"][entry]["category_aro_name"] not in resist_mech[value[hsp]["model_name"]]:
                                                    resist_mech[value[hsp]["model_name"]].append(value[hsp]["ARO_category"][entry]["category_aro_name"])

                                # Drug classes classification
                                        elif value[hsp]["ARO_category"][entry]["category_aro_class_name"] == "Drug Class":
                                            dc += 1
                                            if value[hsp]["model_name"] not in drug_class:
                                                drug_class[value[hsp]["model_name"]] = [value[hsp]["ARO_category"][entry]["category_aro_name"]]
                                            else:
                                                if value[hsp]["ARO_category"][entry]["category_aro_name"] not in drug_class[value[hsp]["model_name"]]:
                                                    drug_class[value[hsp]["model_name"]].append(value[hsp]["ARO_category"][entry]["category_aro_name"])

                                # Gene Family classification
                                        elif value[hsp]["ARO_category"][entry]["category_aro_class_name"] == "AMR Gene Family":
                                            gf += 1
                                            if value[hsp]["model_name"] not in gene_family:
                                                gene_family[value[hsp]["model_name"]] = [value[hsp]["ARO_category"][entry]["category_aro_name"]]
                                            else:
                                                if value[hsp]["ARO_category"][entry]["category_aro_name"] not in gene_family[value[hsp]["model_name"]]:
                                                    gene_family[value[hsp]["model_name"]].append(value[hsp]["ARO_category"][entry]["category_aro_name"])

                                    # Flag to exclude model if it doesn't have classification for rm, gf, or dc
                                    if any(x == 0 for x in [rm, gf, dc]):
                                        del tophits[topmodel]
                                        if topmodel not in excluded:
                                            excluded.append(topmodel)
                                        try:
                                            del resist_mech[topmodel]
                                        except:
                                            pass
                                        try:
                                            del gene_family[topmodel]
                                        except:
                                            pass
                                        try:
                                            del drug_class[topmodel]
                                        except:
                                            pass
                                        # print(jsonfile, hsp)
                                        # print("NOTE: %s excluded because it is missing complete categorization information." % (topmodel))

                            except Exception as e:
                                print(e)
                        else:
                            print("Loose hit encountered. Not being added.")
                # Populates the matrix of typehits
                genes[accession] = tophits
                for x in tophits:
                    if x not in genelist:
                        genelist.append(x)
            # except Exception as e:
            #     print(e)
            except:
                pass

        for e in excluded:
            print("NOTE: %s excluded because it is missing complete categorization information." % (e))

        genelist = sorted(genelist)

        # Create a dictionary that will convert type of hit to num. value
        conversion = {"Perfect": 2, "Strict": 1}

        # Apply conversion so hit criteria is number based
        for sample in genes:
            for gene in genes[sample]:
                genes[sample][gene] = conversion[genes[sample][gene]]
            for thing in genelist:
                if thing not in genes[sample]:
                    genes[sample][thing] = 0

        # Create dataframe from genes dictionary
        if not genes:
            logger.error("Error: No data recovered from JSONs, cannot build heatmap. "
            "Please check you are using RGI results from ver 4.0.0 or greater.")
            exit()
        df = pd.DataFrame.from_dict(genes)

        # Fixed colourmap values (purple, teal, yellow)
        cmap_values = [0, 1, 2, 3]
        custom_cmap = matplotlib.colors.ListedColormap(['#4c0057', '#00948f', '#feed00'])
        norm = matplotlib.colors.BoundaryNorm(cmap_values, custom_cmap.N)

        # If the classification option chosen:
        if self.classification:
            global ax0
            if self.classification == "drug_class":
                df = self.create_categories(drug_class, df)
            elif self.classification == "resistance_mechanism":
                df = self.create_categories(resist_mech, df)
            elif self.classification == "gene_family":
                df = self.create_categories(gene_family, df)

            # Create 3 series, one for each classification type
            class_df1 = self.create_class_series(drug_class, "drug_class")
            class_df2 = self.create_class_series(resist_mech, "resistance_mechanism")
            class_df3 = self.create_class_series(gene_family, "gene_family")

            # Combine the 3 Series into a dataframe with all classification info
            complete_class_df = pd.concat([class_df1, class_df2, class_df3], axis=1, sort=True)

            # Function if possible
            if self.classification == "drug_class":
                classification = "Drug Class"
                complete_class_df= complete_class_df.set_index(["resistance_mechanism", "gene_family"], append=True)["drug_class"].apply(pd.Series).stack()
                complete_class_df= complete_class_df.reset_index()
                complete_class_df.columns = ["model_name", "resistance_mechanism", "gene_family", "number", "drug_class"]
                complete_class_df= complete_class_df.set_index("model_name")
                complete_class_df= complete_class_df.drop(["number"], axis=1)
            elif self.classification == "resistance_mechanism":
                classification = "Resistance Mechanism"
                complete_class_df= complete_class_df.set_index(["drug_class", "gene_family"], append=True)["resistance_mechanism"].apply(pd.Series).stack()
                complete_class_df= complete_class_df.reset_index()
                complete_class_df.columns = ["model_name", "drug_class", "gene_family", "number", "resistance_mechanism"]
                complete_class_df= complete_class_df.set_index("model_name")
                complete_class_df= complete_class_df.drop(["number"], axis=1)
            elif self.classification == "gene_family":
                classification = "AMR Gene Family"
                complete_class_df= complete_class_df.set_index(["drug_class", "resistance_mechanism"], append=True)["gene_family"].apply(pd.Series).stack()
                complete_class_df= complete_class_df.reset_index()
                complete_class_df.columns = ["model_name", "drug_class", "resistane_mechanism", "number", "gene_family"]
                complete_class_df= complete_class_df.set_index("model_name")
                complete_class_df= complete_class_df.drop(["number"], axis=1)

            # Create unique identifiers again for the classifications dataframe
            new_index = []
            counted = {}
            for i,v in enumerate(list(complete_class_df.index.values)):
                if v in counted:
                    counted[v] += 1
                    new_index.append(v+"_"+str(counted[v]))
                else:
                    counted[v] = 0
                    new_index.append(v+"_0")

            # Assign new column to dataframe called uID with unique identifiers
            complete_class_df = complete_class_df.assign(uID=new_index)
            complete_class_df = complete_class_df.reset_index().set_index("uID")
            complete_class_df = complete_class_df.sort_values(by=[self.classification, 'model_name'])
            s = complete_class_df.loc[:,self.classification]
            unique_ids = list(complete_class_df.index.values)
            df = df.reindex(index=unique_ids)

            # Figure parameters if frequency option chosen
            if self.frequency:
                df, freq_dict = self.create_frequency_df(df, self.output)
                df = df.reindex(index=unique_ids)
                df_for_merging = df.copy()
                # Write csv before removing unique IDs
                if not self.cluster:
                    file_name = '{}-{}'.format(self.output, str(len(jsons)))
                    self.write_csv(self.classification, df, s, file_name)
                # Set index as labels, not unique IDs
                df = df.set_index("index")
                column_order = list(df)

                if self.cluster == "samples":
                    df = self.cluster_data(self.cluster, df)
                    # df = df.set_index("index", append=True)
                    # Remove unique ids, cluster, add back unique ids, merge
                    unique_id_values = df_for_merging.index.values
                    df_for_merging = df_for_merging.set_index("index")
                    df_for_merging = df_for_merging.iloc[:,clustered_col]
                    test_df = df_for_merging.assign(TEMP_COLUMN = unique_id_values)
                    test_df = test_df.reset_index().set_index("TEMP_COLUMN")
                    file_name = '{}-{}'.format(self.output, str(len(jsons)))
                    self.write_csv(self.classification, test_df, s, file_name)
                    column_order = list(df)
                elif self.cluster == "both" or self.cluster == "genes":
                    logger.error("Error: Unable to cluster genes because the categorization option was chosen. No heatmap will be generated. Closing program now.")
                    exit()

                # Set the figure size
                fig_width,fig_length,fig,figsize = self.get_figure_dimensions(jsons, unique_ids)

                # Try to draw plot with default sizing
                ax0,ax1,ax2,gs = self.create_plot('cf', 4)

                # Adjust the dimensions
                while True:
                    if self.get_axis_size(fig,ax0)[1] > 150:
                        fig_length = fig_length/2
                        figsize = (fig_width, fig_length)
                        desired_width = (self.get_axis_size(fig,ax0)[1])/3
                        figsize = (desired_width, fig_length)
                        fig = plt.figure(figsize = figsize)
                        ax0,ax1,ax2,gs = self.create_plot('cf', 4)
                    if self.get_axis_size(fig,ax0)[0] < (self.get_axis_size(fig,ax0)[1])/3:
                        # print('LESS THAN 1/3')
                        fig_length = fig_length/2
                        figsize = (fig_width, fig_length)
                        desired_width = (self.get_axis_size(fig,ax0)[1])/3
                        figsize = (desired_width, fig_length)
                        fig = plt.figure(figsize = figsize)
                        ax0,ax1,gs = self.create_plot('c', 4)
                    if self.get_axis_size(fig,ax0)[0] < 10:
                        # print('BASE AXIS TOO SMALL')
                        try:
                            desired_width = desired_width*2
                        except:
                            desired_width = fig_width
                        figsize = (desired_width, fig_length)
                        fig = plt.figure(figsize = figsize)
                        ax0,ax1,gs = self.create_plot('c', 4)
                    if self.get_axis_size(fig,ax0)[1] < 150:
                        if self.get_axis_size(fig,ax0)[0] > (self.get_axis_size(fig,ax0)[1])/3:
                            if self.get_axis_size(fig,ax0)[0] > 10:
                                break

                # Calculate correct categories dimensions to use
                ratio_to_use = (self.get_axis_size(fig,ax0)[0])/8

                if figsize[1] > 100:
                    sns.set(font_scale=1.7)
                if df.shape[0] > 200:
                    sns.set(font_scale=1.0)
                sns.set_style("white")

                ax0,ax1,ax2,gs = self.create_plot('cf', ratio_to_use)

                # Create the heatmap
                g = sns.heatmap(df, cmap=custom_cmap, cbar=False, ax=ax0, norm=norm) #linewidth=0.5
                plt.setp(g.yaxis.get_ticklabels(), rotation=0, fontsize='xx-large')
                plt.setp(g.xaxis.get_ticklabels(), visible=False)
                g.tick_params(bottom=False)
                g.set_xlabel(" ")
                g.yaxis.set_label_position("left")
                g.set_ylabel(" ")
                plt.setp(ax1.get_yticklabels(), visible=False)
                plt.setp(ax1.get_xticklabels(), visible=False)
                ax1_dim = self.get_axis_size(fig, ax1)
                cat, ranges = self.calculate_categories(s, ax1_dim[0])
                # Draw categories
                self.draw_categories(ax1, ranges, cat, ax0, self.display, df)

                # Draw barplot
                self.draw_barplot(freq_dict,ax2, column_order)

                # Save figure
                gs.tight_layout(fig)
                print("Rendering EPS")
                plt.savefig(file_name + '.eps', bbox_inches="tight", format="eps", pad_inches=0.5)
                print("Rendering PNG")
                plt.savefig(file_name + '.png', bbox_inches="tight", format="png", pad_inches=0.5)
                if self.cluster == "samples":
                    print('Output file {fn}: AMR genes categorised by {c} and only unique '
                    'resistome profiles are displayed with ther frequency and have been '
                    'clustered hierarchically (see SciPy documentation). Yellow '
                    'represents a perfect hit, teal represents a strict hit, purple '
                    'represents no hit. Genes with asterisks (*) appear multiple times '
                    'because they belong to more than one {c} category in the '
                    'antibiotic resistance ontology (ARO).'.format(fn=file_name, c=classification))
                else:
                    print('Output file {fn}: AMR genes categorised by {c} and '
                    'only unique resistome profiles are displayed with ther '
                    'frequency. Yellow represents a perfect hit, teal represents '
                    'a strict hit, purple represents no hit. Genes with asterisks '
                    '(*) appear multiple times because they belong to more than '
                    'one {c} category in the antibiotic resistance ontology (ARO).'
                    .format(fn=file_name, c=classification))

            # Categories, but no frequency
            else:
                # Modifies dataframe if cluster option chosen
                if self.cluster == "samples":
                    df_copy = df.drop(["index"], axis=1)
                    df_copy = self.cluster_data(self.cluster, df_copy)
                    df = df.set_index("index", append=True)
                    df = df.iloc[:, clustered_col]
                    df = df.reset_index().set_index("uID")
                elif self.cluster == "both" or self.cluster == "genes":
                    logger.error("Error: Unable to cluster genes because the categorization option was chosen. No heatmap will be generated. Closing program now.")
                    exit()

                # Write csv before removing unique IDs
                file_name = '{}-{}'.format(self.output, str(len(jsons)))
                self.write_csv(self.classification, df, s, file_name)
                df = df.set_index("index")

                # Set the dimension parameters
                fig_width,fig_length,fig,figsize = self.get_figure_dimensions(jsons, unique_ids)

                # Try to draw plot with default sizing
                ax0,ax1,gs = self.create_plot('c', 4)

                # Adjust the dimensions
                while True:
                    if self.get_axis_size(fig,ax0)[1] > 150:
                        fig_length = fig_length/2
                        figsize = (fig_width, fig_length)
                        desired_width = (self.get_axis_size(fig,ax0)[1])/3
                        figsize = (desired_width, fig_length)
                        fig = plt.figure(figsize = figsize)
                        ax0,ax1,gs = self.create_plot('c', 4)
                    if self.get_axis_size(fig,ax0)[0] < (self.get_axis_size(fig,ax0)[1])/3:
                        # print('LESS THAN 1/3')
                        fig_length = fig_length/2
                        figsize = (fig_width, fig_length)
                        desired_width = (self.get_axis_size(fig,ax0)[1])/3
                        figsize = (desired_width, fig_length)
                        fig = plt.figure(figsize = figsize)
                        ax0,ax1,gs = self.create_plot('c', 4)
                    if self.get_axis_size(fig,ax0)[0] < 10:
                        try:
                            desired_width = desired_width*2
                        except:
                            desired_width = fig_width
                        figsize = (desired_width, fig_length)
                        fig = plt.figure(figsize = figsize)
                        ax0,ax1,gs = self.create_plot('c', 4)
                    if self.get_axis_size(fig,ax0)[1] < 150:
                        if self.get_axis_size(fig,ax0)[0] > (self.get_axis_size(fig,ax0)[1])/3:
                            if self.get_axis_size(fig,ax0)[0] > 10:
                                break

                # Calculate correct categories dimensions to use
                ratio_to_use = (self.get_axis_size(fig,ax0)[0])/8

                if figsize[1] > 100:
                    sns.set(font_scale=1.7)
                if df.shape[0] > 500:
                    sns.set(font_scale=0.5)
                elif df.shape[0] > 200:
                    sns.set(font_scale=1.0)
                sns.set_style("white")
                ax0,ax1,gs = self.create_plot('c', ratio_to_use)

                # Create the heatmap
                g = sns.heatmap(df, cmap=custom_cmap, cbar=False, ax=ax0, norm=norm) #linewidth=0.5
                plt.setp(g.yaxis.get_ticklabels(), rotation=0, fontsize='xx-large')
                plt.setp(g.xaxis.get_ticklabels(), rotation=90, fontsize='xx-large')
                plt.setp(ax1.get_yticklabels(), visible=False)
                plt.setp(ax1.get_xticklabels(), visible=False)
                g.set_ylabel(" ")
                g.set_xlabel(" ")
                ax1_dim = self.get_axis_size(fig, ax1)

                cat, ranges = self.calculate_categories(s, ax1_dim[0])
                self.draw_categories(ax1, ranges, cat, ax0, self.display, df)

                # Save figure
                gs.tight_layout(fig)
                print("Rendering EPS")
                plt.savefig(file_name + '.eps', bbox_inches="tight", format="eps", pad_inches=0.5)
                print("Rendering PNG")
                plt.savefig(file_name + '.png', bbox_inches="tight", format="png", pad_inches=0.5)
                if self.cluster == "samples":
                    print('Output file {n}: AMR genes categorised by {c} and '
                    'samples have been clustered hierarchically (see SciPy '
                    'documentation). Yellow represents a perfect hit, teal '
                    'represents a strict hit, purple represents no hit. Genes '
                    'with asterisks (*) appear multiple times because they belong '
                    'to more than one {c} category in the antibiotic resistance '
                    'ontology (ARO).'.format(n=file_name, c=classification))
                else:
                    print('Output file %s: AMR genes categorised by %s. Yellow '
                    'represents a perfect hit, teal represents a strict hit, purple '
                    'represents no hit. Genes with asterisks (*) appear multiple times '
                    'because they belong to more than one %s category in the '
                    'antibiotic resistance ontology (ARO).' %(file_name, classification, classification))

        # No categories
        else:
            if self.frequency:
                df,freq_dict = self.create_frequency_df(df, self.output)
                column_order = list(df)

                if self.cluster:
                    df = self.cluster_data(self.cluster, df)
                    column_order = list(df)

                file_name = '{}-{}'.format(self.output, str(len(jsons)))
                # Write matrix to csv
                df.index.name='gene'
                df.to_csv('{}.csv'.format(file_name))

                # Set the dimension parameters
                fig_width,fig_length,fig,figsize = self.get_figure_dimensions(jsons, genelist)

                # Try to draw plot with default sizing
                if figsize[1] > 100:
                    sns.set(font_scale=1.7)
                # if df.shape[0] > 200:
                #     sns.set(font_scale=1.0)
                sns.set_style("white")
                ax0,ax2,gs = self.create_plot('f', 0)

                # Create the heatmap
                g = sns.heatmap(df, cmap=custom_cmap, cbar=False, ax=ax0, norm=norm) #linewidth=0.5
                plt.setp(g.yaxis.get_ticklabels(), rotation=0, fontsize='xx-large')
                plt.setp(g.xaxis.get_ticklabels(), visible=False)
                g.tick_params(bottom=False)
                g.yaxis.set_label_position("left")
                g.set_ylabel(" ")
                g.set_xlabel(" ")

                # Draw barplot
                self.draw_barplot(freq_dict,ax2, column_order)

                # Save figure
                try:
                    gs.tight_layout(fig)
                except:
                    # print("fixing")
                    # Increase width of plot to avoid matplotlib tight_layout error
                    new_figsize = (fig_width*3, fig_length)
                    fig = plt.figure(figsize = new_figsize)
                    # Try to draw plot with default sizing
                    ax0,ax2,gs = self.create_plot('f', 0)
                    temp_w, temp_h = self.get_axis_size(fig, ax0)
                    # Check if width too small
                    if temp_h/temp_w > 7:
                        new_figsize = (fig_width*5, fig_length)
                        fig = plt.figure(figsize = new_figsize)
                        ax0,ax2,gs = self.create_plot('f', 0)
                    if figsize[1] > 100:
                        sns.set(font_scale=1.7)
                    # if df.shape[0] > 200:
                    #     sns.set(font_scale=1.0)
                    sns.set_style("white")
                    df,freq_dict = self.create_frequency_df(df, self.output)

                    # Create the heatmap
                    g = sns.heatmap(df, cmap=custom_cmap, cbar=False, ax=ax0, norm=norm) #linewidth=0.5
                    plt.setp(g.yaxis.get_ticklabels(), rotation=0, fontsize='xx-large')
                    plt.setp(g.xaxis.get_ticklabels(), visible=False)
                    g.tick_params(bottom=False)
                    g.yaxis.set_label_position("left")
                    g.set_ylabel(" ")
                    g.set_xlabel(" ")

                    # Draw barplot
                    self.draw_barplot(freq_dict,ax2, column_order)
                    gs.tight_layout(fig)

                print("Rendering EPS")
                plt.savefig(file_name + '.eps', bbox_inches="tight", format="eps", pad_inches=0.5)
                print("Rendering PNG")
                plt.savefig(file_name + '.png', bbox_inches="tight", pad_inches=0.5, format="png")
                if self.cluster == 'samples':
                    print('Output file %s: AMR genes are listed in alphabetical order and unique '
                    'resistome profiles are displayed with their frequency have been '
                    'clustered hierarchically (see SciPy documentation). Yellow '
                    'represents a perfect hit, teal represents a strict hit, purple '
                    'represents no hit.' %(file_name))
                elif self.cluster == 'genes':
                    print('Output file %s: AMR genes have been clustered hierarchically '
                    '(see SciPy documentation) and unique resistome profiles are '
                    'displayed with ther frequency. Yellow represents a perfect hit, teal represents a strict hit, purple '
                    'represents no hit.' %(file_name))
                elif self.cluster == 'both':
                    print('Output file %s: AMR genes and unique resistome profiles '
                    'displayed with ther frequency have been clustered hierarchically '
                    '(see SciPy documentation). Yellow represents a perfect hit, teal represents a strict hit, purple '
                    'represents no hit.' %(file_name))
                else:
                    print('Output file %s: AMR genes are listed in alphabetical order and unique '
                    'resistome profiles are displayed with their frequency. Yellow '
                    'represents a perfect hit, teal represents a strict hit, purple '
                    'represents no hit.' %(file_name))

            else:
                # No categories or frequency
                if self.cluster:
                    df = self.cluster_data(self.cluster, df)

                file_name = '{}-{}'.format(self.output, str(len(jsons)))

                # Write matrix to csv
                df.index.name='gene'
                df.to_csv('{}.csv'.format(file_name))

                # Set the dimension parameters
                fig_width,fig_length,fig,figsize = self.get_figure_dimensions(jsons, genelist)

                if figsize[1] > 100:
                    sns.set(font_scale=1.7)
                if fig_width > 200 and fig_length > 200:
                    figsize = (fig_width/2, fig_length/2)
                    fig = plt.figure(figsize = figsize)
                    sns.set(font_scale=1.2)

                sns.set_style("white")

                # Create the heatmap
                g = sns.heatmap(df, cmap=custom_cmap, cbar=False, norm=norm) #linewidth=0.5
                plt.setp(g.yaxis.get_ticklabels(), rotation=0, fontsize='xx-large')
                plt.setp(g.xaxis.get_ticklabels(), rotation=90, fontsize='xx-large')
                g.set_ylabel(" ")
                g.set_xlabel(" ")

                # Save figure
                print("Rendering EPS")
                plt.savefig(file_name + '.eps', bbox_inches="tight", format="eps", pad_inches=0.5)
                print("Rendering PNG")
                plt.savefig(file_name + '.png', bbox_inches="tight", format="png", pad_inches=0.5)
                if self.cluster == 'samples':
                    print('Output file %s: AMR genes are listed in alphabetical order '
                    'and samples have been clustered hierarchically (see SciPy documentation). '
                    'Yellow represents a perfect hit, teal represents a strict hit, purple '
                    'represents no hit.' %(file_name))
                elif self.cluster == 'genes':
                    print('Output file %s: AMR genes have been clustered hierarchically. '
                    'Yellow represents a perfect hit, teal represents a strict hit, purple '
                    'represents no hit.' %(file_name))
                elif self.cluster == 'both':
                    print('Output file %s: AMR genes and samples have been clustered hierarchically '
                    '(see SciPy documentation). Yellow represents a perfect hit, teal represents a strict hit, purple '
                    'represents no hit.' %(file_name))
                else:
                    print('Output file %s: Yellow represents a perfect hit, '
                    'teal represents a strict hit, purple represents no hit.' %(file_name))
