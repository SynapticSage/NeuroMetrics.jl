using GoalFetchAnalysis
using Timeshift
using Timeshift.checkpoint
using Table

# Posthoc adding labels to namped-tuples
@time I, S, F = load_mains(),
                load_shuffles(),
                load_fields();

@time Id = Table.to_dataframe(I, key_name=["shift"], explode=false);
@time Sd = Table.to_dataframe(S,       key_name=["shift"], explode=false);
@time Fd = Table.to_dataframe(F,       key_name=["shift"], explode=false);

save_mains(Id,    overwrite=true)
save_shuffles(Sd, overwrite=true)
save_fields(Fd,   overwrite=true)
