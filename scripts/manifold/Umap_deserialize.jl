using DataStructures: OrderedDict

using GoalFetchAnalysis
import Utils.namedtup: ntopt_string
using  DataFramesMeta
using UMAP
using Serialization


embedding, embedding2, inds, animal, day = deserialize(datadir("manifold","ca1pfc_manifolds.serial"))

@time spikes, beh, ripples, cells = Load.load(animal, day);
Rca1, Rpfc = (Munge.spiking.torate(@subset(spikes,:area .== ar), beh) for ar in ("CA1","PFC"))
nsamp = min(100_000, size(Rca1,1))
