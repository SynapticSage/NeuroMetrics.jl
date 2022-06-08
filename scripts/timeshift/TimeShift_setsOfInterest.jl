# Highest priority
marginals_highprior = [
    ["x", "y"],
    ["currentHeadEgoAngle", "currentPathLength"],
    ["currentHeadEgoAngle", "currentPathLength", "stopWell"],
    ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
]


marginals_superhighprior_shuffle = [
    ["x", "y"],
]

# Priority for shuffle
# Remove the super high dimensional 5D marginal -- would take 50 days to compute
marginals_highprior_shuffle = [
    ["x", "y"],
    ["currentHeadEgoAngle", "currentPathLength"],
    ["currentHeadEgoAngle", "currentPathLength", "stopWell"],
]

# Lower Priority
marginals_lowprior = [
    ["currentHeadEgoAngle"],
    ["currentPathLength"],
    ["stopWell"],
    ["x", "y", "currentHeadEgoAngle"],
    ["x", "y", "stopWell"],
    ["x", "y", "currentHeadEgoAngle", "currentPathLength"]
]

# Lower Priority
marginals_lowprior_shuffle = [
    ["currentHeadEgoAngle"],
    ["currentPathLength"],
    ["stopWell"],
    ["x", "y", "currentHeadEgoAngle"],
    ["x", "y", "stopWell"],
]
