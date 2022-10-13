
using Munge.dynamic
examples = get_groupedexamples(spikes, beh)
templates = get_templates(examples)
dtwtab = get_dtwtable(examples, templates)
Xwarp = apply_warps(X, dtwtab)
Twarp = warped_df_to_tensor(Xwarp, [:startWell, :stopWell], :data;
                            inner_dim_names=:neuron);
