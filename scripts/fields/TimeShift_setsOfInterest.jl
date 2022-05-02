# Highest priority
marginals_highprior = [
    ["x", "y"],
    ["currentHeadEgoAngle", "currentPathLength"],
    ["currentHeadEgoAngle", "currentPathLength", "stopWell"],
    ["x", "y", "currentHeadEgoAngle", "currentPathLength", "stopWell"]
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
