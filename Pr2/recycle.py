    # Minimizing numeric-domains.
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            if start_state[i][j] != "*#":
                if start_state[i][j][0] != "*":
                    for k in range(dimension_of_table):
                        if (str(i) + str(k)) != (str(i) + str(j)):
                            num_domains[str(i) + str(k)].remove(int(start_state[i][j][0]))

                if start_state[i][j][1] != "#":
                    if i == 0:
                        if j == 0:
                            color_domains[str(i + 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j + 1)].remove(start_state[i][j][1])
                        elif j == dimension_of_table - 1:
                            color_domains[str(i) + str(j - 1)].remove(start_state[i][j][1])
                            color_domains[str(i + 1) + str(j)].remove(start_state[i][j][1])
                        else:
                            color_domains[str(i + 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j - 1)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j + 1)].remove(start_state[i][j][1])

                    elif i == dimension_of_table - 1:
                        if j == 0:
                            color_domains[str(i - 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j + 1)].remove(start_state[i][j][1])
                        elif j == dimension_of_table - 1:
                            color_domains[str(i) + str(j - 1)].remove(start_state[i][j][1])
                            color_domains[str(i - 1) + str(j)].remove(start_state[i][j][1])
                        else:
                            color_domains[str(i - 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j - 1)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j + 1)].remove(start_state[i][j][1])

                    else:
                        if j == 0:
                            color_domains[str(i - 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i + 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j + 1)].remove(start_state[i][j][1])
                        elif j == dimension_of_table - 1:
                            color_domains[str(i - 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i + 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j - 1)].remove(start_state[i][j][1])
                        else:
                            color_domains[str(i - 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i + 1) + str(j)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j - 1)].remove(start_state[i][j][1])
                            color_domains[str(i) + str(j + 1)].remove(start_state[i][j][1])

# ----------------------------------------------------------------------------------------------------------------------

from csp import Csp
import copy
import numpy as np


def main():
    start_state, txt_list = read_file()

    dimension_of_table = int(txt_list[0][1])  # vertical or horizontal dimension of table

    # **************** Creating a CSP problem ****************
    # Define numeric constraints of each variable.
    ini_numeric_constraints = numeric_constraints_generator(dimension_of_table)
    # print(ini_numeric_constraints)

    # Define colorful constraints of each variable.
    ini_colorful_constraints = colorful_constraints_generator(dimension_of_table)
    # print(ini_colorful_constraints)

    # Define numeric domain of each variable.
    num_domains = numeric_domain_generator(dimension_of_table)
    # print(num_domains)

    # Define colorful domain of each variable.
    color_domains = color_domain_generator(dimension_of_table, txt_list)
    # print(color_domains)

    num_domains, color_domains, u_vars = limiting_domain(dimension_of_table, start_state, num_domains, color_domains,
                                                         ini_numeric_constraints, ini_colorful_constraints)

    # **************** Generating data for CSP Object ****************

    num_of_constraints = create_num_of_constraints(ini_numeric_constraints, ini_colorful_constraints)

    root_node = Csp(num_of_constraints, ini_colorful_constraints, ini_numeric_constraints,
                    start_state, u_vars, num_domains, color_domains)

    print(num_of_constraints)

    # A dictionary. Keys are variables and values are nested list.
    # First element in each list is number and second is color. (number is int)
    vars_domain = combine_colorful_numeric_domains(num_domains, color_domains)
    print(vars_domain)
    selected_var = mrv(vars_domain, u_vars, num_of_constraints)
    # d = forward_checking(var_num_neighbors, var_color_neighbors, domain_dict, unassigned_vars, [])


'''
Returns a two dimensional array of input text 
that a first index is row and the second is column of input text. 
'''
def read_file():
    file = open('./input.txt')
    text_input = file.read()
    print(text_input)

    file.close()
    txt_list = []
    t = text_input.split("\n")
    for word in t:
        txt_list.append(word.split())  # A list that contains lines of input text.
    start_state = txt_list[2:]  # two dimensional array. (our initial state without two first lines of input txt)
    # print(start_state[0][0][0])
    return start_state, txt_list


# Returns numeric constraints of variables.
def numeric_constraints_generator(dimension_of_table):
    ini_numeric_constraints = {}
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            const_list = []  # list of constraints that a variable participated
            k = 0
            while k < dimension_of_table:
                const_list.append(str(i) + str(k))
                const_list.append(str(k) + str(j))
                k += 1
            # To not add self element to its constraints list
            const_list.remove(str(i) + str(j))
            const_list.remove(str(i) + str(j))
            ini_numeric_constraints.update({str(i) + str(j): copy.deepcopy(const_list)})

    return ini_numeric_constraints


# Returns colorful constraints of variables.
def colorful_constraints_generator(dimension_of_table):
    ini_colorful_constraints = {}
    for i in range(dimension_of_table):

        for j in range(dimension_of_table):
            const_list = []  # list of constraints that a variable participated

            if i == 0:
                if j == 0:
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j + 1))
                elif j == dimension_of_table - 1:
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i + 1) + str(j))
                else:
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i) + str(j + 1))

            elif i == dimension_of_table - 1:
                if j == 0:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i) + str(j + 1))
                elif j == dimension_of_table - 1:
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i - 1) + str(j))
                else:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i) + str(j + 1))
            else:
                if j == 0:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j + 1))
                elif j == dimension_of_table - 1:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                else:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i) + str(j + 1))

            ini_colorful_constraints.update({str(i) + str(j): copy.deepcopy(const_list)})

    return ini_colorful_constraints


# Returns numeric domain of variables.
def numeric_domain_generator(dimension_of_table):
    num_domains = {}
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            num_domains[str(i) + str(j)] = [num for num in range(1, dimension_of_table + 1)]
    return num_domains


# Returns color domain of variables.
def color_domain_generator(dimension_of_table, txt_list):
    color_domains = {}
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            color_domains[str(i) + str(j)] = [color for color in txt_list[1]]
    return color_domains


# Minimizing numeric-domains.
def limiting_domain(dimension_of_table, start_state, num_domains, color_domains, ini_numeric_constraints,
                    ini_colorful_constraints):
    unassigned_variables = []
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            if start_state[i][j] != "*#":
                if (start_state[i][j][0] != "*") and (start_state[i][j][1] != "#"):
                    num_domains.pop(str(i) + str(j))

                if start_state[i][j][0] != "*":

                    if (str(i) + str(j)) in list(num_domains.keys()):
                        num_domains[str(i) + str(j)].clear()
                        num_domains[str(i) + str(j)].append(int(start_state[i][j][0]))

                        for var in ini_numeric_constraints[str(i) + str(j)]:
                            if var in list(num_domains.keys()):
                                if int(start_state[i][j][0]) in num_domains[var]:
                                    num_domains[var].remove(int(start_state[i][j][0]))

                if start_state[i][j][1] != "#":
                    color_domains[str(i) + str(j)].clear()
                    color_domains[str(i) + str(j)].append(start_state[i][j][1])

                    for var in ini_colorful_constraints[str(i) + str(j)]:
                        if start_state[i][j][1] in color_domains[var]:
                            color_domains[var].remove(start_state[i][j][1])

            if start_state[i][j][0] == "*" or start_state[i][j][1] == "#":
                unassigned_variables.append(str(i) + str(j))

    return num_domains, color_domains, unassigned_variables


def create_num_of_constraints(ini_numeric_constraints, ini_colorful_constraints):
    num_of_constraints = {}
    for var in ini_numeric_constraints:
        num_of_constraints[var] = len(ini_numeric_constraints[var])

    for var in ini_colorful_constraints:
        num_of_constraints[var] += len(ini_colorful_constraints[var])
    return num_of_constraints


def combine_colorful_numeric_domains(num_domains, color_domains):
    # for var in num_domains:
    #     num_domains[var] = list(np.array(num_domains[var]).astype(str))
    domain_dict = {}
    somelists = []
    for var in num_domains:
        somelists.append(num_domains[var])
        somelists.append(color_domains[var])
        var_domain = [[a, b] for a in somelists[0] for b in somelists[1]]
        domain_dict.update({var: var_domain})
        somelists.clear()
    return domain_dict


def mrv(domain_dict, unassigned_vars, num_of_constraints):
    tmp_dict = {}
    for var in domain_dict:
        if var in unassigned_vars:
            tmp_dict.update({var: len(domain_dict[var])})

    sorted_dic = {k: v for k, v in sorted(tmp_dict.items(), key=lambda item: item[1])}
    # Key for minimum value of values.
    key_min = min(sorted_dic.keys(), key=(lambda k: sorted_dic[k]))
    # Check if we have two equal vale with different keys or not. If we do, we should use degree heuristic.
    res = sum(x == sorted_dic[key_min] for x in sorted_dic.values())

    # print(tmp_dict)
    if res == 1:
        # print("mrv")
        return key_min
    else:
        # print("degree")
        return degree(num_of_constraints)


def degree(num_of_constraints):
    # Key for minimum value of values.
    key_min = min(num_of_constraints.keys(), key=(lambda k: num_of_constraints[k]))
    return key_min


def forward_checking(var_num_neighbors, var_color_neighbors, domain_dict, unassigned_vars, assigned_val):
    failure = False
    for neighbor in var_num_neighbors:
        if neighbor in unassigned_vars:
            if assigned_val in domain_dict[neighbor]:
                domain_dict[neighbor].remove(assigned_val)

    for neighbor in var_color_neighbors:
        if neighbor in unassigned_vars:
            if assigned_val in domain_dict[neighbor]:
                domain_dict[neighbor].remove(assigned_val)

    for var_domain in domain_dict:
        if len(var_domain) == 0:
            print("mmmm----")
            failure = True
    return failure, domain_dict


if __name__ == "__main__":
    main()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

from csp import Csp
import copy
import numpy as np

COLOR_PRIORITY = {}
DIM = 0


def main():
    start_state, txt_list = read_file()

    dimension_of_table = int(txt_list[0][1])  # vertical or horizontal dimension of table
    global DIM
    DIM = dimension_of_table

    priority = len(txt_list[1])
    for color in txt_list[1]:
        COLOR_PRIORITY.update({color: priority})
        priority -= 1

    # start_state[1][1] = "3g"
    # print(start_state)

    # **************** Creating a CSP problem ****************
    # Define numeric constraints of each variable.
    ini_numeric_constraints = numeric_constraints_generator(dimension_of_table)
    # print(ini_numeric_constraints)

    # Define colorful constraints of each variable.
    ini_colorful_constraints = colorful_constraints_generator(dimension_of_table)
    # print(ini_colorful_constraints)

    # Define numeric domain of each variable.
    num_domains = numeric_domain_generator(dimension_of_table)
    # print(num_domains)

    # Define colorful domain of each variable.
    color_domains = color_domain_generator(dimension_of_table, txt_list)
    # print(color_domains)

    num_domains, color_domains, u_vars = limiting_domain(dimension_of_table, start_state, num_domains, color_domains,
                                                         ini_numeric_constraints, ini_colorful_constraints)

    # print(num_of_constraints)

    # A dictionary. Keys are variables and values are nested list.
    # First element in each list is number and second is color. (number is int)
    vars_domain = combine_colorful_numeric_domains(num_domains, color_domains)
    print(vars_domain)

    num_of_constraints = create_num_of_constraints(ini_numeric_constraints, ini_colorful_constraints)

    f, d = forward_checking(ini_numeric_constraints["12"], ini_colorful_constraints["12"], vars_domain, u_vars, "1g")
    print(d)

    # **************** Generating data for CSP Object ****************

    root_node = Csp(num_of_constraints, ini_colorful_constraints, ini_numeric_constraints,
                    start_state, u_vars, num_domains, color_domains, vars_domain)

    assignment = {}
    for key in vars_domain:
        if key not in u_vars:
            assignment.update({key: start_state[int(key[0])][int(key[1])]})

    backtrack(assignment, (dimension_of_table * dimension_of_table), root_node)

    print(is_consistent("2b", assignment, ini_numeric_constraints["10"], ini_colorful_constraints["10"]))


'''
Returns a two dimensional array of input text 
that a first index is row and the second is column of input text. 
'''


def read_file():
    file = open('./input.txt')
    text_input = file.read()
    print(text_input)

    file.close()
    txt_list = []
    t = text_input.split("\n")
    for word in t:
        txt_list.append(word.split())  # A list that contains lines of input text.
    start_state = txt_list[2:]  # two dimensional array. (our initial state without two first lines of input txt)
    # print(start_state[0][0][0])
    return start_state, txt_list


# Returns numeric constraints of variables.
def numeric_constraints_generator(dimension_of_table):
    ini_numeric_constraints = {}
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            const_list = []  # list of constraints that a variable participated
            k = 0
            while k < dimension_of_table:
                const_list.append(str(i) + str(k))
                const_list.append(str(k) + str(j))
                k += 1
            # To not add self element to its constraints list
            const_list.remove(str(i) + str(j))
            const_list.remove(str(i) + str(j))
            ini_numeric_constraints.update({str(i) + str(j): copy.deepcopy(const_list)})

    return ini_numeric_constraints


# Returns colorful constraints of variables.
def colorful_constraints_generator(dimension_of_table):
    ini_colorful_constraints = {}
    for i in range(dimension_of_table):

        for j in range(dimension_of_table):
            const_list = []  # list of constraints that a variable participated

            if i == 0:
                if j == 0:
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j + 1))
                elif j == dimension_of_table - 1:
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i + 1) + str(j))
                else:
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i) + str(j + 1))

            elif i == dimension_of_table - 1:
                if j == 0:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i) + str(j + 1))
                elif j == dimension_of_table - 1:
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i - 1) + str(j))
                else:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i) + str(j + 1))
            else:
                if j == 0:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j + 1))
                elif j == dimension_of_table - 1:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                else:
                    const_list.append(str(i - 1) + str(j))
                    const_list.append(str(i + 1) + str(j))
                    const_list.append(str(i) + str(j - 1))
                    const_list.append(str(i) + str(j + 1))

            ini_colorful_constraints.update({str(i) + str(j): copy.deepcopy(const_list)})

    return ini_colorful_constraints


# Returns numeric domain of variables.
def numeric_domain_generator(dimension_of_table):
    num_domains = {}
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            num_domains[str(i) + str(j)] = [num for num in range(1, dimension_of_table + 1)]
    return num_domains


# Returns color domain of variables.
def color_domain_generator(dimension_of_table, txt_list):
    color_domains = {}
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            color_domains[str(i) + str(j)] = [color for color in txt_list[1]]
    return color_domains


# Minimizing numeric-domains.
def limiting_domain(dimension_of_table, start_state, num_domains, color_domains, ini_numeric_constraints,
                    ini_colorful_constraints):
    unassigned_variables = []
    for i in range(dimension_of_table):
        for j in range(dimension_of_table):
            if start_state[i][j] != "*#":

                if start_state[i][j][0] != "*":
                    num_domains[str(i) + str(j)].clear()
                    num_domains[str(i) + str(j)].append(int(start_state[i][j][0]))

                    for var in ini_numeric_constraints[str(i) + str(j)]:
                        if int(start_state[i][j][0]) in num_domains[var]:
                            num_domains[var].remove(int(start_state[i][j][0]))

                if start_state[i][j][1] != "#":
                    color_domains[str(i) + str(j)].clear()
                    color_domains[str(i) + str(j)].append(start_state[i][j][1])

                    for var in ini_colorful_constraints[str(i) + str(j)]:
                        if start_state[i][j][1] in color_domains[var]:
                            color_domains[var].remove(start_state[i][j][1])

            if start_state[i][j][0] == "*" or start_state[i][j][1] == "#":
                unassigned_variables.append(str(i) + str(j))

    return num_domains, color_domains, unassigned_variables


def create_num_of_constraints(ini_numeric_constraints, ini_colorful_constraints):
    num_of_constraints = {}
    for var in ini_numeric_constraints:
        num_of_constraints[var] = len(ini_numeric_constraints[var])

    for var in ini_colorful_constraints:
        num_of_constraints[var] += len(ini_colorful_constraints[var])
    return num_of_constraints


def combine_colorful_numeric_domains(num_domains, color_domains):
    # for var in num_domains:
    #     num_domains[var] = list(np.array(num_domains[var]).astype(str))
    domain_dict = {}
    somelists = []
    for var in num_domains:
        somelists.append(num_domains[var])
        somelists.append(color_domains[var])
        var_domain = [(str(a) + str(b)) for a in somelists[0] for b in somelists[1]]
        domain_dict.update({var: var_domain})
        somelists.clear()
    return domain_dict


def mrv(domain_dict, unassigned_vars, num_of_constraints):
    tmp_dict = {}
    for var in domain_dict:
        if var in unassigned_vars:
            tmp_dict.update({var: len(domain_dict[var])})

    sorted_dic = {k: v for k, v in sorted(tmp_dict.items(), key=lambda item: item[1])}
    # Key for minimum value of values.
    key_min = min(sorted_dic.keys(), key=(lambda k: sorted_dic[k]))
    # Check if we have two equal vale with different keys or not. If we do, we should use degree heuristic.
    res = sum(x == sorted_dic[key_min] for x in sorted_dic.values())

    # print(tmp_dict)
    if res == 1:
        # print("mrv")
        return key_min
    else:
        # print("degree")
        return degree(num_of_constraints)


def degree(num_of_constraints):
    # Key for minimum value of values.
    key_min = min(num_of_constraints.keys(), key=(lambda k: num_of_constraints[k]))
    return key_min


def forward_checking(var_num_neighbors, var_color_neighbors, domain_dict, unassigned_vars, assigned_val):
    failure = False
    # print(var_color_neighbors)
    # print(unassigned_vars)
    # print(assigned_val)
    # if assigned_val in
    # print(domain_dict["00"])

    # Release limits on numbers.
    for neighbor in var_num_neighbors:
        if neighbor in unassigned_vars:
            tmp_list = [subl for subl in domain_dict[neighbor] if subl[0] != assigned_val[0]]
            domain_dict[neighbor] = copy.deepcopy(tmp_list)

            # # Priority color limitation by numeric dependency.
            # if int(assigned_val[0]) == DIM:
            #     t_list = [value for value in domain_dict[neighbor]
            #               if COLOR_PRIORITY[value[1]] < COLOR_PRIORITY[assigned_val[1]]]
            #     domain_dict[neighbor] = copy.deepcopy(t_list)
            #
            # elif int(assigned_val[0]) == 1:
            #     t_list = [value for value in domain_dict[neighbor]
            #               if COLOR_PRIORITY[value[1]] > COLOR_PRIORITY[assigned_val[1]]]
            #     domain_dict[neighbor] = copy.deepcopy(t_list)

    # Release limits on colors.
    for neighbor in var_color_neighbors:
        if neighbor in unassigned_vars:
            tmp_list = [subl for subl in domain_dict[neighbor] if subl[1] != assigned_val[1]]
            domain_dict[neighbor] = copy.deepcopy(tmp_list)

            # Priority color limitation colorful dependency.
            if COLOR_PRIORITY[assigned_val[1]] == max(COLOR_PRIORITY.keys(), key=(lambda k: COLOR_PRIORITY[k])):
                t_list = [value for value in domain_dict[neighbor] if int(value[0]) < int(assigned_val[0])]
                domain_dict[neighbor] = copy.deepcopy(t_list)
            elif COLOR_PRIORITY[assigned_val[1]] == min(COLOR_PRIORITY.keys(), key=(lambda k: COLOR_PRIORITY[k])):
                t_list = [value for value in domain_dict[neighbor] if int(value[0]) > int(assigned_val[0])]
                domain_dict[neighbor] = copy.deepcopy(t_list)

            less_list = [value for value in domain_dict[neighbor] if (int(value[0]) < int(assigned_val[0]))
                         and (COLOR_PRIORITY[value[1]] < COLOR_PRIORITY[assigned_val[1]])]
            more_list = [value for value in domain_dict[neighbor] if (int(value[0]) > int(assigned_val[0]))
                         and (COLOR_PRIORITY[value[1]] > COLOR_PRIORITY[assigned_val[1]])]
            f_list = less_list + more_list
            domain_dict[neighbor] = copy.deepcopy(f_list)

    for key in domain_dict:
        if len(domain_dict[key]) == 0:
            failure = True
    return failure, domain_dict


def backtrack(assignment, dimension_of_table, csp):
    if len(assignment) == dimension_of_table:
        return assignment
    variable = mrv(csp.domains_dict, csp.unassigned_vars, csp.num_of_constraints)
    for value in csp.domains_dict[variable]:
        if is_consistent(value, assignment, csp.numeric_const_dict[variable], csp.colorful_const_dict[variable]):
            pass


def is_consistent(value, assignment, var_num_neighbors, var_color_neighbors):
    for var in assignment:
        if (var in var_num_neighbors) and (assignment[var][0] == value[0]):
            print("aaaaa")
            return False
        if (var in var_color_neighbors) and (assignment[var][1] == value[1]):
            print("jjjjjj")
            return False
        if var in var_color_neighbors:
            if assignment[var][0] != "*":
                if COLOR_PRIORITY[assignment[var][1]] > COLOR_PRIORITY[value[1]] \
                        and (int(assignment[var][0]) < int(value[0])):
                    print("ccccc")
                    return False
        if var in var_color_neighbors:
            if assignment[var][0] != "*":
                if COLOR_PRIORITY[assignment[var][1]] < COLOR_PRIORITY[value[1]] \
                        and (int(assignment[var][0]) > int(value[0])):
                    print("ddddd")
                    return False

    return True


if __name__ == "__main__":
    main()
# //////////////////////////////////////////////////////////////////////////////////////////////////////////// -2/5/2021

def backtrack(assignment_list, dimension_of_table, csp):
    print(csp.state, "KKKK")
    if len(assignment_list) == dimension_of_table:
        return assignment_list
    variable = mrv(csp.domains_dict, csp.unassigned_vars, csp.num_of_constraints)
    for value in csp.domains_dict[variable]:
        if is_consistent(value, assignment_list, csp.state, csp.numeric_const_dict[variable],
                         csp.colorful_const_dict[variable]):

            assignment_list.update({variable: value})
            domain_copied = copy.deepcopy(csp.domains_dict)
            is_empty, new_domain = forward_checking(csp.numeric_const_dict[variable], csp.colorful_const_dict[variable],
                                                    domain_copied, csp.unassigned_vars, value)
            if not is_empty:
                v = variable
                ns = copy.deepcopy(csp.state)  # New state
                ns[int(variable[0])][int(variable[1])] = value
                child_node = Csp(csp.num_of_constraints.pop(v), csp.colorful_const_dict.pop(v),
                                 csp.numeric_const_dict.pop(v), ns, csp.unassigned_vars.remove(variable),
                                 csp.domain_dict_num.pop(v), csp.domain_dict_color.pop(v), new_domain.pop(v))

                result = backtrack(assignment_list, dimension_of_table, child_node)
                if result != "failure":
                    return result
            assignment_list.pop(variable)
            return "failure"