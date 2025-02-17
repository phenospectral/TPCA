import numpy as np
from dash.exceptions import PreventUpdate
from scipy import stats
import plotly.express as px
from statsmodels.regression.linear_model import OLS
from xhtml2pdf import pisa
import pandas as pd
from dash import Dash, dcc, html, Input, Output, callback, dash_table, ctx
from dash.dash_table.Format import Format, Scheme, Trim


class TPCA:

    """This class provides functions to import proteomic and transcriptomic gene-product differential expression data,
    join them on the basis of GeneID, and analyse statistical relationships between transcript and protein responses.
    Hierachical functional catgeory labels of the form TopLevel.SecondLevel.ThirdLevel.FourthLevel... and so on,
    provided in a column of the joined data may be used to selectively plots certain gene bins"
    """

    def __init__(self):
        
        self.app = Dash()
        self.transcript_data_file = 'data/Peng_Setaria_Protein_Data_20240521.xlsx'
        self.protein_data_file = 'data/Peng_Setaria_Transcript_Data_20240521.xlsx'
        self.protein_data = pd.DataFrame()
        self.transcript_data = pd.DataFrame()
        self.protein_transcript_data = pd.DataFrame()
        self.level_fbins = {}
        self.graphs = []
        self.output_df= pd.DataFrame()
        self.fbin_transcript_protein_correlation_table = {}
        self.fbin_gene_tree_data={
            "index": [],
            "name": [],
            "parent": [],
            "value": [],
            "level":[],
            "r": [],
            "r-squared": [],
            "nobs": [],
        }
        self.config = config = {

            # provide the path to the proteomic data table (*.xlsx)
            'protein_data_file': 'data/Peng_Setaria_Protein_Data_20240521.xlsx',

            # provide the path to the transcriptomic data table (*.xlsx)
            'transcript_data_file': 'data/Peng_Setaria_Transcript_Data_20240521.xlsx',

            # name of the unique geneID column in the proteomic data (for joining to the transcript data via exact match)
            # in proteomic data, this will be the specific or (in multi-protein protein groups),
            # the geneID of the most well-evidenced (razor) protein
            'protein_data_index_col': 'GeneID',

            # name of the unique geneID column in the transcript data (for joining to the protein data via exact match)
            'transcript_data_index_col': 'GeneID',

            # name of the column in the protein data containing the list of gene IDs represented in the protein group
            'protein_data_gene_ids_col': 'GeneIDs',

            # name of the column in the protein data containing the protein abundance values to be used in visualisations
            'protein_data_abundance_col': 'Intensity',

            # name of the column in the protein data containing the log-transformed protein response data
            # (e.g. log2(FoldChange) values
            'protein_data_logFC_response_data_col': 'log2FC_protein',

            # name of the column in the transcript data containing the log-transformed transcript response data
            # (e.g. log2(FoldChange) values
            'transcript_data_logFC_response_data_col': 'log2FC_transcript',

            # name of the column in the transcript data containing the log-transformed transcript response data
            # (e.g. log2(FoldChange) values
            'transcript_data_FC_data_col': 'foldChange',

            # name of the column in the protein data table that contains the protein group expression ratio column
            'protein_data_FC_data_col': 'Ratio H/L',

            # name of the column in the protein data containing the full annotation / description of the columns
            'protein_data_description_col': 'Fasta headers',

            # name of the column (in either transcript or protein data table) that contains the hierarchical functional bin labels
            'fbin_col': 'Functional Bin',

            # name of the column (in either transcript or protein data table) that contains the hierarchical functional bin labels
            'gene_id_col': 'newlocusName',

            # name of the column in the protein data table that contains the peptide count column
            'protein_data_peptide_count_col': 'Peptides',

            # name of the column in the protein data table that contains the protein group Score column
            'protein_data_score_col': 'Score',

            # name of the column in the protein data table that contains the protein group Score column
            'transcript_data_basemean_col': 'baseMean',

            # specify the figure height in pixels
            'figure_height': 800

        }

    def get_parent_fbin(self, fbin):

        if fbin.lower() == "all":
            parent_fbin = ""
            return parent_fbin

        else:
            subfbins = fbin.split(".")
            # join the first n subfbins together again to make a string that represents the current level depth into the hierarchy
            parent_fbin = ".".join(subfbins[0:len(subfbins)-1])
            if parent_fbin == "":
                parent_fbin = "all"

        return parent_fbin

    def configure(self, config):

        for key in config.keys():
            if key in self.config.keys():
                self.config[key] = config[key]

    def get_level_fbins_from_protein_transcript_data(self):
        """
            from the protein_transcript dataframe, pull out the non-redundant
            list of functional bins in the column with header "NAME"
            """
        fbins = ["all"]
        fbins.extend(self.protein_transcript_data[self.config["fbin_col"]].unique())

        """
        prepare a dictionary that will store lists of functionally categories 
        at each level in the hierarchy of functional categories and return it
        """

        self.level_fbins = {}

        # iterate through each of the unique functional category bins (fbins)
        for fbin in fbins:

            # split the fbin into its level components (subfbins) by splitting wherever there is a . character
            subfbins = fbin.split(".")

            # for each hierarchical level represented by the current fbin
            for level in range(0, len(subfbins)+1):

                # join the first n subfbins together again to make a string that represents the current level depth into the hierarchy
                level_subfbin = ".".join(subfbins[0:level])
                if level_subfbin != "":

                    # if we don't yet have a key in the level_bins dictionary under which to store the list of fbins in that level, create an
                    # empty array under a key representing the current level
                    if level not in self.level_fbins.keys():
                        self.level_fbins[level] = []

                    # if we haven't already added the current level_subfbin under its level key (we may have), do so now.
                    if level_subfbin not in self.level_fbins[level]:
                        self.level_fbins[level].append(level_subfbin)

                        self.fbin_gene_tree_data["index"].append(level_subfbin)
                        self.fbin_gene_tree_data["name"].append(level_subfbin)
                        self.fbin_gene_tree_data["parent"].append(self.get_parent_fbin(level_subfbin))
                        self.fbin_gene_tree_data["value"].append(0)
                        self.fbin_gene_tree_data["level"].append(level)
                        self.fbin_gene_tree_data["r"].append(0)
                        self.fbin_gene_tree_data["r-squared"].append(0)
                        self.fbin_gene_tree_data["nobs"].append(0)

    def get_fbin_gene_tree_data_from_protein_transcript_data(self):

        for index, row in self.protein_transcript_data.iterrows():
            parent_fbin = row[self.config["fbin_col"]]
            if parent_fbin != "":
                level = len(parent_fbin.split("."))+2
                self.fbin_gene_tree_data["index"].append(row[self.config["protein_data_description_col"]])
                self.fbin_gene_tree_data["name"].append(row[self.config["protein_data_description_col"]])
                self.fbin_gene_tree_data["parent"].append(parent_fbin)
                self.fbin_gene_tree_data["value"].append(row[self.config["protein_data_abundance_col"]])
                self.fbin_gene_tree_data["level"].append(level)
                self.fbin_gene_tree_data["r"].append(0)
                self.fbin_gene_tree_data["r-squared"].append(0)
                self.fbin_gene_tree_data["nobs"].append(1)

    def prepare(self):

        self.get_protein_data()
        self.get_transcript_data()

        # Map (join) the transcript data rows onto the ends of their respective protein rows on the basis of index_col match
        self.protein_transcript_data = self.protein_data.join(self.transcript_data, how='inner')

        # Output the joined protein-transcript data frame to an Excel file
        self.protein_transcript_data.to_excel("data/Protein-Transcript_Data.xlsx")

        # Map the joined protein-transcript dataset rows out into a dict of arrays under unique level_fbin keys
        self.get_level_fbins_from_protein_transcript_data()
        self.get_fbin_gene_tree_data_from_protein_transcript_data()
        self.fbin_gene_tree_data = pd.DataFrame(self.fbin_gene_tree_data).astype({'value': 'float64', 'r': 'float64', 'r-squared': 'float64', 'level': 'int64', 'nobs': 'int64'})
        self.fbin_gene_tree_data.set_index("index", inplace=True)

    def report_block_template(self, report_type, graph_url, caption=''):

        graph_block = ""

        if report_type == 'interactive':
            graph_block = '<iframe style="border: none;" src="exports/charts/{graph_url}.embed" width="100%" height="600px"></iframe>'

        elif report_type == 'static':
            graph_block = (''
                '<a href="{graph_url}" target="_blank">' # Open the interactive graph when you click on the image
                    '<img style="height: 400px;" src="exports/charts/{graph_url}.png">'
                '</a>')

        report_block = ('' + graph_block +
                '{caption}' + # Optional caption to include below the graph
                '<br>'      + # Line break
                '<a href="{graph_url}" style="color: rgb(190,190,190); text-decoration: none; font-weight: 200;" target="_blank">'+
                    'Click to comment and see the interactive graph' + # Direct readers to Plotly for commenting, interactive graph
                '</a>' +
                '<br>' +
                '<hr>') # horizontal line

        return report_block.format(graph_url=graph_url, caption=caption)

    # function to carry out ordinary least squares analysis on two arrays of values and return the result
    def ols(self, x, y):
        model = OLS(y, x)
        results = model.fit()
        return results

    def analyse_subset_association(self, subset):

        r = stats.pearsonr(subset[self.config["transcript_data_logFC_response_data_col"]], subset[self.config["protein_data_logFC_response_data_col"]])[0]
        rsquared = r ** 2
        r_formatted = "{:.2f}".format(r)
        rsquared_formatted = "{:.2f}".format(rsquared)
        nobs = len(subset)

        return {'r': r, 'r_formatted': r_formatted, 'r-squared': rsquared, 'r-squared_formatted': rsquared_formatted, 'nobs': nobs}

    def get_subset_by_fbin(self, fbin):
        if fbin.lower() == "all":
            subset = self.protein_transcript_data
        else:
            print(fbin)
            subset = self.protein_transcript_data[self.protein_transcript_data[self.config["fbin_col"]].str.contains(fbin)]
        return subset

    def plot_scatter(self, functional_bin, write_html):

        """
        function to generate an interactive scatterplot of transcript responses against protein responses for a
        subset of genes/proteins belonging to a particular functional category at any level. It generates a .html
        file and returns ordinary least squares regression analysis summary results to the caller.

        It takes two positional arguments:

        'protein_transcript_data' (pandas DataFrame) = the protein:transcript data pandas dataframe returned by joining the
        transcript and protein dataframes on geneID match.

        'functional_bin' (str) = the name of the functional bin. If the value "All" is passed, the function will generate
        a scatter plot of all genes/proteins in the set instead of just those in a specific functional bin.

        """

        if functional_bin is None:
            functional_bin = "all"

        styles = {
            'pre': {
                'border': 'thin lightgrey solid',
                'overflowX': 'scroll'
            }
        }

        subset = self.get_subset_by_fbin(functional_bin)

        if len(subset) > 1:
            results = self.analyse_subset_association(subset)
            if results is not None:

                title = ("log2FC(Transcript) vs log2FC(Protein): " +
                         "[" + functional_bin + "]" +
                         " (r = " + str(round(results['r'], 3)) +
                         ", r2 = " + str(round(results['r-squared'], 3)) +
                         ")")

                scatter_fig = px.scatter(subset,
                                         x=self.config["transcript_data_logFC_response_data_col"],
                                         y=self.config["protein_data_logFC_response_data_col"],
                                         color=self.config["protein_data_abundance_col"],
                                         size=self.config["protein_data_abundance_col"],
                                         trendline='ols',
                                         color_continuous_scale=px.colors.sequential.Bluered,
                                         title=title,
                                         hover_name=self.config["protein_data_description_col"],
                                         custom_data=[self.config["gene_id_col"]],
                                         template='plotly_white')

                scatter_fig.update_layout(
                    plot_bgcolor='white',
                    hovermode="x unified",
                    hoverlabel=dict(bgcolor='rgba(180, 180, 180, 0.5)', font_size=10, font_color='rgba(0, 0, 0, 1)'),
                    clickmode='event+select',
                    height=self.config["figure_height"]
                )

                if write_html:
                    html_file = functional_bin.replace(".", "_").replace("/", "-").replace("\\", "-") + ".html"
                    html_file_fullpath = "exports/charts/" + html_file

                    scatter_fig.write_html(html_file_fullpath)

                return scatter_fig

            else:
                return None
        else:
            return None

    # Utility function
    def convert_html_to_pdf(self, source_html, output_filename):
        # open output file for writing (truncated binary)
        result_file = open(output_filename, "w+b")

        try:
            # convert HTML to PDF
            pisa_status = pisa.CreatePDF(
                source_html,  # the HTML to convert
                dest=result_file)  # file handle to receive result

            # close output file
            result_file.close()  # close output file

            # # return True on success and False on errors
            if pisa_status.err == False:
                return result_file
        except:
            print("An error occurred")

    def get_protein_data(self):

        # read the file
        self.protein_data = pd.read_excel(self.config['protein_data_file'], engine="calamine")
        self.protein_data.set_index(self.config["protein_data_index_col"], inplace=True)

        # Clean up the proteome data by removing junk rows:
        self.protein_data.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.protein_data = self.protein_data[~self.protein_data[self.config["protein_data_gene_ids_col"]].str.contains("REV__")]
        self.protein_data = self.protein_data[~self.protein_data[self.config["protein_data_gene_ids_col"]].str.contains("CON__")]
        self.protein_data = self.protein_data[self.protein_data[self.config["protein_data_peptide_count_col"]]>1]
        self.protein_data = self.protein_data[self.protein_data[self.config["protein_data_FC_data_col"]].notnull()]
        self.protein_data = self.protein_data[self.protein_data[self.config["protein_data_score_col"]]>10]
        self.protein_data = self.protein_data[self.protein_data[self.config["protein_data_score_col"]].notnull()]

        # Output the cleaned data to an Excel file:
        self.protein_data.to_excel('data/Protein_Data_Cleaned.xlsx')

    def get_transcript_data(self):

        #read the file
        self.transcript_data = pd.read_excel(self.config['transcript_data_file'], engine="calamine")
        self.transcript_data.set_index(self.config["transcript_data_index_col"], inplace=True)

        # Clean up the transcript data by removing junk rows:
        self.transcript_data.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.transcript_data = self.transcript_data[self.transcript_data[self.config["transcript_data_FC_data_col"]].notnull()]
        self.transcript_data = self.transcript_data[self.transcript_data[self.config["transcript_data_basemean_col"]]>10]
        self.transcript_data = self.transcript_data[self.transcript_data[self.config["transcript_data_logFC_response_data_col"]].notnull()]

        # Output the cleaned data to an Excel file:
        self.transcript_data.to_excel('data/Transcript_Data_Cleaned.xlsx')

    def analyse_TP_fbins(self, plotted_level_bins, options):

        # loop through all levels of the fbin hierarchy
        for level in range(1, len(self.level_fbins)):

            # loop through all the fbins at the current level
            for fbin in self.level_fbins[level]:

                if "all" in plotted_level_bins or fbin in plotted_level_bins:

                    # run the function to generate a scatter plot for the current fbin and collect the returned results from the OLS analysis

                    subset = self.get_subset_by_fbin(fbin)
                    results = None

                    if len(subset)>2:
                        results = self.analyse_subset_association(subset)
                        if options['write_html']:
                            self.plot_scatter(fbin, True)

                    # if OLS results are returned, store them in the dictionary of named arrays
                    if results is not None:

                        self.fbin_gene_tree_data.at[fbin, 'level'] = level
                        self.fbin_gene_tree_data.at[fbin, 'r'] = results["r"]
                        self.fbin_gene_tree_data.at[fbin, 'r-squared'] = results["r-squared"]
                        self.fbin_gene_tree_data.at[fbin, 'nobs'] = results["nobs"]

        self.fbin_gene_tree_data.to_excel("exports/charts/Transcript_Protein_Corr.xlsx")

        return

    def plot(self, options):

        styles = {
            'pre': {
                'border': 'thin lightgrey solid',
                'overflowX': 'scroll'
            }
        }

        self.app = Dash("Transcripto-Proteomic Correlation Analyzer")

        sunburst_fig = px.sunburst(self.fbin_gene_tree_data,
                                   names="name",
                                   parents="parent",
                                   values="value",
                                   color="r",
                                   color_continuous_scale="RdBu",
                                   color_continuous_midpoint=0)

        sunburst_fig.update_layout(height=self.config["figure_height"])

        title = ("log2FC(Transcript) vs log2FC(Protein) [All Functional Bins]: " +
                 " (r = " + str(round(self.fbin_gene_tree_data.at["all", "r"], 3)) +
                 " r2 = " + str(round(self.fbin_gene_tree_data.at["all", "r-squared"], 3)) +
                 ")")

        scatter_fig = px.scatter(self.protein_transcript_data,
                                 x=self.config["transcript_data_logFC_response_data_col"],
                                 y=self.config["protein_data_logFC_response_data_col"],
                                 color=self.config["protein_data_abundance_col"],
                                 size=self.config["protein_data_abundance_col"],
                                 trendline='ols',
                                 color_continuous_scale=px.colors.sequential.Bluered,
                                 title=title,
                                 hover_name=self.config["protein_data_description_col"],
                                 custom_data=[self.config["gene_id_col"]],
                                 template='plotly_white')

        scatter_fig.update_layout(
            plot_bgcolor='white',
            hovermode="x unified",
            hoverlabel=dict(bgcolor='rgba(180, 180, 180, 0.5)', font_size=10, font_color='rgba(0, 0, 0, 1)'),
            clickmode='event+select',
            height=self.config["figure_height"]
        )

        fbin_gene_tree_data_cols = [{"name": i, "id": i} for i in self.fbin_gene_tree_data.columns]

        for c in range(len(fbin_gene_tree_data_cols)):
            if fbin_gene_tree_data_cols[c]['name'] in ['r', 'r-squared']:
                fbin_gene_tree_data_cols[c]['type'] = 'numeric'
                fbin_gene_tree_data_cols[c]['format'] = Format(precision=3, scheme=Scheme.decimal)

        self.app.layout = html.Div([

            html.Div([
                html.Div([
                    dcc.Graph(id='sunburst', figure=sunburst_fig),


                ], style={'padding': 10, 'flex': 1}),

                html.Div([
                    dcc.Graph(id='scatter', figure=scatter_fig)

                ], style={'padding': 10, 'flex': 1})
            ], style={'display': 'flex', 'flexDirection': 'row'}),

            dcc.Tabs([

                dcc.Tab(label='Functional Bins', children=[

                    dash_table.DataTable(self.fbin_gene_tree_data.to_dict('records'),
                                         columns=fbin_gene_tree_data_cols,
                                         id='fbin-table',
                                         sort_action="native",
                                         page_size=20,
                                         )
                    ]),
                dcc.Tab(label='Gene Products', children=[

                    dash_table.DataTable(self.protein_transcript_data.to_dict('records'),
                                         columns=[{"name": i, "id": i} for i in self.protein_transcript_data.columns if
                                                  i in [self.config["protein_data_gene_ids_col"],
                                                        self.config["protein_data_description_col"]]],
                                         id='gene-table',
                                         sort_action="native",
                                         page_size=20,
                                         ),
                    ]),
                # dcc.Tab(label='Gene Info', children=[
                #
                # ]),
            ]),
        ])

        self.callbacks(self.app)
        self.app.run(debug=True)


    def callbacks(self, app):

        @app.callback(
            Output('scatter', 'figure'),
            Input('sunburst', 'clickData'),
            prevent_initial_call = True
        )

        def update_scatter_on_sunburst_click(clickData):
            if clickData is not None:
                clicked_point = clickData["points"][0]
                if clicked_point["label"] == clicked_point["entry"]:
                    if clicked_point["label"] != "all":
                        target_fbin = clicked_point["parent"]
                    else:
                        target_fbin = "all"
                elif clicked_point["label"] == "":
                    target_fbin = "all"
                else:
                    target_fbin = clicked_point["label"]

                new_figure = self.plot_scatter(target_fbin, False)
                if new_figure is not None:
                    return new_figure
                else:
                    raise PreventUpdate

        @app.callback(
            Output('gene-table', 'data'),
            Input('scatter', 'selectedData'),
            Input('scatter', 'clickData'),

            prevent_initial_call=True,
            suppress_callback_exceptions=True,
        )
        def display_selected_gene_products(selectedData, clickData):

            gene_set = []

            triggered_prop_id = ctx.triggered[0]["prop_id"]
            if triggered_prop_id == "scatter.selectedData":
                for point in selectedData["points"]:
                    gene_set.append(point["customdata"][0])

            elif triggered_prop_id == "scatter.clickData":
                gene_set = []
                gene_set.append(clickData['points'][1]["customdata"][0])

            selected_subset_data = self.protein_transcript_data[self.protein_transcript_data[self.config["gene_id_col"]].isin(gene_set)].to_dict("records")

            return selected_subset_data

        @app.callback(
            Output('fbin-table', 'data'),
            Input('scatter', 'selectedData'),
            Input('scatter', 'clickData'),
            Input('sunburst', 'clickData'),

            prevent_initial_call=True,
            suppress_callback_exceptions=True,
        )
        def display_selected_fbins(selectedData, scatterclickdata, sunburstclickdata):

            gene_set = []

            triggered_prop_id = ctx.triggered[0]["prop_id"]

            if ("scatter." in triggered_prop_id):
                if triggered_prop_id == "scatter.selectedData":
                    for point in selectedData["points"]:
                        gene_set.append(point["customdata"][0])

                elif triggered_prop_id == "scatter.clickData":
                    gene_set = []
                    gene_set.append(scatterclickdata['points'][1]["customdata"][0])

                selected_subset_pt_data = self.protein_transcript_data[self.protein_transcript_data[self.config["gene_id_col"]].isin(gene_set)]
                fbins = selected_subset_pt_data[self.config["fbin_col"]].unique()
                selected_fbin_gene_tree_data = self.fbin_gene_tree_data[self.fbin_gene_tree_data['name'].isin(fbins)].to_dict("records")
                return selected_fbin_gene_tree_data

            elif triggered_prop_id == "sunburst.clickData":
                if sunburstclickdata is not None:
                    clicked_point = sunburstclickdata["points"][0]
                    if clicked_point["label"] == clicked_point["entry"]:
                        if clicked_point["label"] != "all":
                            target_fbin = clicked_point["parent"]
                        else:
                            target_fbin = "all"
                    elif clicked_point["label"] == "":
                        target_fbin = "all"
                    else:
                        target_fbin = clicked_point["label"]

                    if target_fbin == "all":
                        selected_fbin_gene_tree_data = self.fbin_gene_tree_data.to_dict("records")
                    else:
                        selected_fbin_gene_tree_data = self.fbin_gene_tree_data[self.fbin_gene_tree_data['name'].str.contains(target_fbin)].to_dict("records")

                    return selected_fbin_gene_tree_data