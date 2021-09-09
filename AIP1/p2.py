from Child_node_bfs import Child
import copy
import time

DEPTH = 0


def main():
    # a = [1, "kk"]
    # c = Child("ali", 5, a)
    # print(c.parent)
    # print(c.action)
    # print(c.state)
    start = time.time()
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

    iterative_deepening_search(start_state, number_of_columns)

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


def recursive_dls(node, limit, ncolumns, frontier, explored):

    if goal_test(node.state, ncolumns):
        solution(node, len(explored), len(explored) + len(frontier))
        # print(node.state)
        # print(len(explored) + len(frontier))
        return "True"
    elif limit == 0:
        # print(node.state, node.n, limit)
        # print("------------------------")
        return "CutOff"
    else:
        explored.append(node)
        frontier.pop(0)
        cutoff_occeured = False
        # print(len(frontier))
        # print(node.actions)
        # print(node.state, node.n, limit)
        # print("------------------------")
        for act in node.actions:
            new_state = rseult(node.state, act, node.actions)
            child = Child(new_state, node.state, node, actions(new_state), act)
            frontier.append(child.state)
            result = recursive_dls(child, (limit-1), ncolumns, frontier, explored)
            if result == "CutOff":
                cutoff_occeured = True
            elif result != "failure":
                return result
        if cutoff_occeured:
            return "CutOff"
        else:
            return "failure"


def iterative_deepening_search(start_state, ncol):
    cutoff = 0
    while True:
        res = depth_limited_search(start_state, cutoff, ncol)
        if res == "True":
            return res
        elif res == "CutOff":
            cutoff += 1
            # print(cutoff, "t")
            # print(cutoff)
        elif res == "failure":
            print("Failed")


def depth_limited_search(initial_state, limit, ncol):

    actions_of_state = actions(initial_state)
    root_node = Child(initial_state, "", "", actions_of_state, "")
    frontier = [root_node.state, ]
    explored = []
    rs = recursive_dls(root_node, limit, ncol, frontier, explored)
    return rs


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
            print("Action done to this state: ",
                  "move card {}".format(ans1_node.state[int(action[0])][len(ans1_node.state[int(action[0])]) - 1]),
                  "from row{}".format(int(action[0])+1), "to row{}".format(int(action[1])+1))
            print("----------------------------------------------------------")
        else:
            print("Final State:", ans1_node.state)
            print("----------------------------------------------------------")

        # solution_dic.update({str(ans1_node.state): action})

    print("Number of generated nodes:{}".format(generated_nodes))
    print("Number of explored nodes:{}".format(expanded_nodes))


if __name__ == "__main__":
    main()
