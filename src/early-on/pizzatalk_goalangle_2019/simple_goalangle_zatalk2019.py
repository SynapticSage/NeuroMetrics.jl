# IPython log file

X=pd.read_csv('RY9_69_04_CBi.1DeepCut_resnet101_goalmaze_implantAug31shuffle1_100000filtered.csv', header=[1,2], index_col=[0])

X = X.swaplevel(0,1, axis=1)
X.likelihood.hist()
likely_points = (X.likelihood > 0.1).all(axis=1)
print(len(X))
X = X[likely_points]
print(len(X))

X = X.swaplevel(0,1, axis=1)
sns.jointplot(X['anterior_drive'].x, X['anterior_drive'].y, kind='scatter')
likely_points.mean()
y
x  = X['left_drive'] + (X['left_drive'] - X['right_drive'])/2
x.name='midpoint'


X = pd.read_csv('RY9_69_04_CBi.1DeepCut_resnet101_goalmaze_implantAug31shuffle1_100000filtered.csv')
X
X=pd.read_csv('RY9_69_04_CBi.1DeepCut_resnet101_goalmaze_implantAug31shuffle1_100000filtered.csv', header=[1,2], index_col=[0])
X
X = X.swaplevel(0,1, axis=1)
X
X.likelihood.hist()
X.likelihood > 0.1
likely_points = (X.likelihood > 0.1).all(axis=1)
get_ipython().run_line_magic('magic', '')
get_ipython().run_line_magic('mkdir', '/home/ryoung/Projects/behavior/goalangle')
get_ipython().run_line_magic('logon', '/home/ryoung/Projects/behavior/goalangle/simple_goalangle_zatalk2019.py')
get_ipython().run_line_magic('logstart/home/ryoung/Projects/behavior/goalangle/simple_goalangle_zatalk2019.py', '')
get_ipython().run_line_magic('logstart', '/home/ryoung/Projects/behavior/goalangle/simple_goalangle_zatalk2019.py')
likely_points = (X.likelihood > 0.1).all(axis=1)
X = X[likely_points.index]
X
X.index.to_frame()
likely_points
X = X[likely_points]
X
len(x)
len(X)
sns.jointplot(X['anterior_drive'].x, X['anterior_drive'].y, kind='scatter')
X = X.swaplevel(0,1, axis=1)
sns.jointplot(X['anterior_drive'].x, X['anterior_drive'].y, kind='scatter')
likely_points.mean()
ZX
X
X['mid-point'] = X['left_drive'] - X['right_drive']
X['left-drive']
X
X['mid_point'] = X['left_drive'] + (X['left_drive'] - X['right_drive'])/2
X.loc[:, 'mid_point'] = X['left_drive'] + (X['left_drive'] - X['right_drive'])/2
X
X['left_drive'] + (X['left_drive'] - X['right_drive'])/2
X = X.swaplevel(0,1, axis=1)
X.drop(columns='likelihood', inplace=True)
X
X.loc[0:1, 'mid_point'] = X['left_drive'] + (X['left_drive'] - X['right_drive'])/2
XX
X = X.swaplevel(0,1, axis=1)
X.loc[0:1, 'mid_point'] = X['left_drive'] + (X['left_drive'] - X['right_drive'])/2
X
x, y  = X['left_drive'] + (X['left_drive'] - X['right_drive'])/2
x
y
x  = X['left_drive'] + (X['left_drive'] - X['right_drive'])/2
x
x.x
x.name='midpoint'
x
X = X.assign(x)
X = X.concat(x,axis=1)
X = pd.concat((X,x), axis=1)
X
get_ipython().run_line_magic('edit', '/home/ryoung/Projects/behavior/goalangle/simple_goalangle_zatalk2019.py')
get_ipython().run_line_magic('edit', '/home/ryoung/Projects/behavior/goalangle/simple_goalangle_zatalk2019.py')
X
x
get_ipython().run_line_magic('edit', '/home/ryoung/Projects/behavior/goalangle/simple_goalangle_zatalk2019.py')
x.columns
x.columns.to_list()
get_ipython().run_line_magic('logon', '')
get_ipython().run_line_magic('edit', '/home/ryoung/Projects/behavior/goalangle/simple_goalangle_zatalk2019.py')
get_ipython().run_line_magic('logoff', '')
logong
get_ipython().run_line_magic('logon', '')
new_index = [('midpoint', x) for x in x.columns.to_list()]
x.columns = pd.MultiIndex.from_tuples(new_index)
x
pd.concat((X,x), axis=1)
X = pd.concat((X,x), axis=1)
get_ipython().run_line_magic('logstate', '')
get_ipython().run_line_magic('edit', '/home/ryoung/Projects/behavior/goalangle/simple_goalangle_zatalk2019.py')
X
get_ipython().run_line_magic('logoff', '')
get_ipython().run_line_magic('logstop', '')
get_ipython().run_line_magic('pinfo', 'logstart')
appen()
appen()
get_ipython().run_line_magic('pinfo', 'logstart')
app()
get_ipython().run_line_magic('logstart', '/home/ryoung/Projects/behavior/goalangle/simple_goalangle_zatalk2019.py append')
new_index = [('midpoint', x) for x in x.columns.to_list()]
x.columns = pd.MultiIndex.from_tuples(new_index)
X = pd.concat((X,x), axis=1)
X
goalvector = X.anterior_drive - X.midpoint
goalvector.name = 'goalvector'
goalvector.columns = pd.MultiIndex.from_tuples([('goalvector',x) for x in goalvector.columns])
goalvector
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('who', '')
goalvector
goalvector.name='heading-vector'
goalvector.columns
goalvector = goalvector.rename(columns={'goalvector','heading'},axis=1,level=0)
goalvector = goalvector.rename(columns={'goalvector','heading'},level=0)
goalvector = goalvector.rename(columns={'goalvector':'heading'},level=0)
goalvector
X = pd.concat((X,goalvector), axis=1)
get_ipython().run_line_magic('logstate', '')
X.heading.hist()
plt.suptitle('Heading distributions')
pd.DataFrame([[1,2],[3,4]], columns=[('a'),('b')], index=np.arange(2))
get_ipython().run_line_magic('logoff', '')
testa= pd.DataFrame([[1,2],[3,4]], columns=[('a'),('b')], index=np.arange(2))
testb=  pd.DataFrame([[1,2],[3,4]], columns=[('a',1),('a',2)], index=np.arange(2))
goal_locations = pd.DataFrame([
[605, 117], [1004, 140], [1026, 542], [612, 548]],
columns=['x','y'],
index=pd.Index([1,2,3,4], name='goal_identity'))
goal_locations
for gl, data in goal_locations.iterrows():
    print(gl, data)
    
for gl, data in goal_locations.iterrows():
    print(gl, data)
    
    
goal_vectors = []
for gl, data in goal_locations.iterrows():
    print(gl, data)
    
    goal_heading = X.heading[['x','y']] - data
    goal_heading.name = gl
    
goal_vectors = []
for gl, data in goal_locations.iterrows():
    print(gl, data)
    
    goal_heading = X.heading[['x','y']] - data
    goal_heading.name = gl
    goal_vectors.append(goal_heading)
    
goal_vectors
get_ipython().run_line_magic('logoff', '')
