using GoalFetchAnalysis
using Timeshift
using Timeshift.checkpoint
using Table

# Posthoc adding labels to namped-tuples
@time I, S, F = load_mains(),
          load_shuffles(),
          load_fields();

Id = Table.to_dataframe(I, keyname=["shift"])
Sd = to_dataframe(S)
Fd = to_dataframe(F)

save_mains(Id, overwrite=true)
save_shuffles(Sd, overwrite=true)
save_fields(Fd, overwrite=true)
