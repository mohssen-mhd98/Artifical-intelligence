from Child_node import Child
import copy
import time

DEPTH = 0


def main():
    start = time.time()

    file = open('./input.txt')
    text_input = file.read()
    print(text_input)
    print("----------------------------------------------------------")
    file.close()
    txt_list = []
    t = text_input.split("\n")
    number_of_columns = int(t[0][4])
    number_of_rows = int(t[0][0])
    for word in t:
        txt_list.append(word.split())
    start_state = txt_list[1:]
    # print(heuristic(start_state, number_of_columns))
    a_star(start_state, number_of_columns)

    stop = time.time()
    print("Time to solve\n", round(stop - start, 3), " seconds")


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
        pass
    else:
        return new_state


def goal_test(state, ncolumns):
    flag = False

    for row in state:
        if len(row) != ncolumns:
            if row[0] != "#":
                return flag
    for row in state:
        if row[0] != "#":
            color = row[0][1]
            for element in row:
                if element[1] != color:
                    return flag

    for row in state:
        if row[0] != "#":
            for idx in range(1, len(row)):
                if int(row[idx-1][0]) < int(row[idx][0]):
                    return flag
    return not flag


def a_star(initial_state, ncolumns):

    global EXPANDED_NODE
    actions_of_state = actions(initial_state)
    root_node = Child(initial_state, "", "", actions_of_state, "", heuristic(initial_state, ncolumns))
    path_cost = 0
    # print(actions_of_state)
    if goal_test(initial_state, ncolumns):
        solution(root_node, 0, 1)
        return True
    frontier = [[heuristic(root_node.state, ncolumns), root_node], ]  # {f(node): node}
    explored = []
    frontier_state = [[heuristic(root_node.state, ncolumns), root_node.state], ]

    while True:
        if len(frontier) == 0:
            print("No Answer!")
            return False

        node = frontier.pop(0)[1]
        fs = frontier_state.pop(0)[1]
        if node.state != fs:
            print("state is not equal to node state!")
            break

        if goal_test(node.state, ncolumns):
            solution(node, len(explored), len(explored) + len(frontier_state))
            return True
        explored.append(node.state)
        actions_of_state = actions(node.state)
        # print(node.state)
        for act in actions_of_state:
            flag = True
            new_state = rseult(node.state, act, actions_of_state)
            # print(new_state)
            child = Child(new_state, node.state, node, actions(new_state), act, heuristic(new_state, ncolumns))  # path_cost child = path_cost parent + 1
            if child.state not in explored:  # or child.state not in frontier_state
                if child.state not in frontier_state:
                    # print(new_state)
                    frontier.append([child.path_cost + heuristic(child.state, ncolumns), child])
                    frontier_state.append([child.path_cost + heuristic(child.state, ncolumns), child.state])
                    # print(child.state, frontier[frontier_state.index([child.path_cost + heuristic(child.state, ncolumns), child.state])][0])

        frontier_state.sort(key=lambda x: int(x[0]))
        frontier.sort(key=lambda x: int(x[0]))


def heuristic(state, ncolumns):
    num_of_diff_colors = 0
    distant = 0
    different_color = []
    # print(state)
    r = 0
    for row in state:
        if ncolumns > len(row) > 1:
            r += 1
        for col in row:
            if col[0] == "#":
                continue
            if col[1] not in different_color:
                different_color.append(col[1])
            exact_point = abs(ncolumns - int(col[0]))
            curr_point = row.index(col)
            if curr_point != exact_point:
                distant += 1  # abs(ncolumns - curr_point) + exact_point
        if row[0][0] != "#":
            tmp = len(different_color) - 1
            num_of_diff_colors += tmp
            different_color.clear()
        # print(row, "\n", distant, num_of_diff_colors, different_color)
    total_cost = 0.6*distant + 1.4*(num_of_diff_colors + r)  # 0.3 -- 0.7
    # print(state, distant, num_of_diff_colors, r, total_cost/2)
    return total_cost/2


def solution(answer_node, expanded_nodes, generated_nodes):
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
            print("Action done to this state:", str(int(action[0])+1) + " to " + str(int(action[1])+1))
            print("----------------------------------------------------------")
        else:
            print("Final State:", ans1_node.state)
            print("----------------------------------------------------------")

        # solution_dic.update({str(ans1_node.state): action})

    print("Number of generated nodes:{}".format(generated_nodes))
    print("Number of explored nodes:{}".format(expanded_nodes))


if __name__ == "__main__":
    main()
