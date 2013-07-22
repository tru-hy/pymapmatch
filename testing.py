from pymapmatch.route import *

route = Route([[1,1], [1,3], [3,5], [1,2]])

print list(route.hits_in_range((1, 2), 1))
