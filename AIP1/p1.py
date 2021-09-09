from Child_node_bfs import Child
import copy
import time

DEPTH = 0


def main():
    start = time.time()
    # a = [1, "kk"]
    # c = Child("ali", 5, a)
    # print(c.parent)
    # print(c.action)
    # print(c.state)
    file = open('./input.txt')
    text_input = file.read()
    print(text_input)
    print("----------------------------------------------------------")
    file.close()
    txt_list = []
    t = text_input.split("\n")
    number_of_columns = int(t[0][4])
    for word in t:
        txt_list.append(word.split())
    start_state = txt_list[1:]

    breadth_first_search(start_state, number_of_columns)

    stop = time.time()
    print("Time to solve\n", round(stop - start, 3))


def actions(state):
    last_elements = []
    action_list = []
    for s in state:
        if s[-1][0] != "#":
            last_elements.append(int(s[-1][0]))
        else:
            last_elements.append(s[-1][0])

    for i in range(len(last_elements)):
        for idx in range(len(last_elements)):
            if last_elements[idx] != "#" and last_elements[i] != "#":
                if last_elements[idx] > last_elements[i]:
                    action_list.append(str(i) + str(idx))
            elif last_elements[i] != "#":
                action_list.append(str(i) + str(idx))

    return action_list


def rseult(curr_state, action, state_actions):
    new_state = copy.deepcopy(curr_state)

    if action in state_actions:
        ac = action
        src = int(ac[0])
        dst = int(ac[1])

        if new_state[src][-1] != "#":
            el = new_state[src].pop()
            if len(new_state[src]) == 0:
                new_state[src].append("#")
        if new_state[dst][0] == "#":
            new_state[dst].pop()
            new_state[dst].append(el)
        else:
            new_state[dst].append(el)
    if new_state == curr_state:
        print("duplicated state!!")
        pass
    else:
        return new_state


def goal_test(state, ncolumns):
    flag = False

    for row in state:
        if len(row) != ncolumns:
            if row[0] != "#":
                return flag
    # Checks all color are same in one row.
    for row in state:
        if row[0] != "#":
            color = row[0][1]
            for element in row:
                if element[1] != color:
                    return flag
    # Check numbers in a row sorted decreasingly.
    for row in state:
        if row[0] != "#":
            for idx in range(1, len(row)):
                if int(row[idx-1][0]) < int(row[idx][0]):
                    return flag
    return not flag


def breadth_first_search(initial_state, ncolumns):

    global EXPANDED_NODE
    actions_of_state = actions(initial_state)
    root_node = Child(initial_state, "", "", actions_of_state, "")
    path_cost = 0
    # print(actions_of_state)
    if goal_test(initial_state, ncolumns):
        solution(root_node, 0, 1, ncolumns)
        return True
    frontier = [root_node, ]
    explored = []
    frontier_state = [root_node.state, ]

    while True:
        if len(frontier) == 0:
            return False
        node = frontier.pop(0)
        fs = frontier_state.pop(0)
        if node.state != fs:
            print("state is not equal to node state!")
            break
        # print("state:", node.state)
        # print("actions:", node.actions)
        # print("----------------------------")
        explored.append(node.state)
        actions_of_state = actions(node.state)

        for act in actions_of_state:
            new_state = rseult(node.state, act, actions_of_state)
            # print(new_state)
            child = Child(new_state, node.state, node, actions(new_state), act)
            if child.state not in explored:  # or child.state not in frontier_state
                if goal_test(child.state, ncolumns):
                    # print(child.state)
                    # print(len(explored)+len(frontier_state)) # expanded nodes = len(explored)
                    # print(len(explored))
                    solution(child, len(explored), len(explored)+len(frontier_state), ncolumns)
                    return True
                if child.state not in frontier_state:
                    # print(new_state)
                    frontier.append(child)
                    frontier_state.append(child.state)


def solution(answer_node, expanded_nodes, generated_nodes, ncolumns):
    global DEPTH
    DEPTH = 0
    solution_dic = {}
    ans_node = answer_node
    while True:
        if answer_node == "":
            break
        answer_node = answer_node.parent_obj
        DEPTH += 1

    print("Depth of solution is:{}".format(DEPTH-1))

    for i in range(DEPTH):
        ans1_node = ans_node
        for j in range(i+1, DEPTH):
            ans1_node = ans1_node.parent_obj
            if j == DEPTH - 2:
                action = ans1_node.action

        if i == DEPTH - 2:
            action = ans_node.action

        if ans1_node != ans_node:
            print(ans1_node.state)
            print("Action done to this state: ",
                  "move card {}".format(ans1_node.state[int(action[0])][len(ans1_node.state[int(action[0])]) - 1]),
                  "from row{}".format(int(action[0])+1), "to row{}".format(int(action[1])+1))
            print("----------------------------------------------------------")
        else:
            print("Final State:", ans1_node.state)
            print("----------------------------------------------------------")

        # solution_dic.update({str(ans1_node.state): action})

    print("Number of generated nodes:{}".format(generated_nodes))
    print("Number of expanded nodes:{}".format(expanded_nodes))


if __name__ == "__main__":
    main()
