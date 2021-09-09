
class Child:

    def __init__(self, state, parent, parent_obj, actions, action, h):
        self.state = state
        self.parent = parent
        self.parent_obj = parent_obj
        self.actions = actions
        self.action = action
        self.h = h
        if self.parent_obj != "":
            self.path_cost = self.parent_obj.path_cost + 1  # step_cost = 1
        else:
            self.path_cost = 0


# --------------------------------------
# 1) IDS
# class Child:
#
#     def __init__(self, state, parent, parent_obj, actions):
#         self.state = state
#         self.parent = parent
#         self.parent_obj = parent_obj
#         self.actions = actions
#         # self.path_cost = self.parent.path_cost + step_cost()
#         pass

# 2) BFS
# class Child:
#
#     def __init__(self, state, parent, parent_obj, actions, action):
#         self.state = state
#         self.parent = parent
#         self.parent_obj = parent_obj
#         self.actions = actions
#         self.action = action
#         # self.path_cost = self.parent.path_cost + step_cost()
#         pass
