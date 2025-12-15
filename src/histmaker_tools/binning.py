from src.histmaker_tools.utils import neg_format


def bin_quantity(
    df,
    binned_quantity_field,
    binning_field,
    bins,
    output_prefix="binned_E_reco_over_true",
    histogram_labels="Ereco/Etrue;Ereco/Etrue;Events",
):
    """

    :param df: the root df on which to perform the binning
    :param binned_quantity_field: The field in the root dataframe that contains the quantity to be binned (e.g., E_reco/E_true for jet energy resolution)
    :param binning_field: The field in the root dataframe that contains the quantity on which to perform the binning (e.g., True jet energy)
    :param bins: A Python list of bins of the binning_field.
    :param output_prefix: The prefix to use for the binned histograms.
    :param histogram_labels: The labels to use for the histograms for the FCCAnalyses plots.
    :return:
    """
    histograms = []
    bins_formatted = [neg_format(b) for b in bins]
    for i in range(len(bins) - 1):
        df_column_name = output_prefix + "_{}_{}".format(
            bins_formatted[i], bins_formatted[i + 1]
        )
        df = df.Define(
            df_column_name,
            f"FCCAnalyses::Utils::filter_number_by_bin({binned_quantity_field}, {binning_field}, {bins[i]}, {bins[i+1]})",
        )
        current_bin_hist = df.Histo1D(
            (
                df_column_name,
                histogram_labels,
                5000,
                0,
                2.0,
            ),
            df_column_name,
        )
        histograms.append(current_bin_hist)
    return df, histograms
