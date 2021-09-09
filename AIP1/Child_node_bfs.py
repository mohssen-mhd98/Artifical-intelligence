class Child:

    def __init__(self, state, parent, parent_obj, actions, action):
        self.state = state
        self.parent = parent
        self.parent_obj = parent_obj
        self.actions = actions
        self.action = action
        # self.path_cost = self.parent.path_cost + step_cost()