
dic = {'00': 6, '01': 7, '02': 6, '10': 7, '11': 8, '12': 7, '20': 6, '21': 7, '22': 6}
sorted_dic = {k: v for k, v in sorted(dic.items(), key=lambda item: item[1])}
key_min = min(dic.keys(), key=(lambda k: dic[k]))
res = sum(x == dic[key_min] for x in dic.values())
print(sorted_dic)
print(key_min)
print(dic[key_min])
print(res)

somelists = [
   [1, 2, 3],
   ['a', 'b'],
]

cart_prod = [[a, b] for a in somelists[0] for b in somelists[1]]

print(cart_prod)

a = [1, 2]
b = [1, 2]
#
#
# def main():
#     u = [1, 2, 3]
#     i = 2
#     ass(u, i)
#
#
# def ass(u, i):
#     i = i - 1
#     print("up", u)
#     if i == 0:
#         u.pop()
#         print("lll")
#         return True
#     else:
#         ass(u, i)
#     print("down", u)
#
#
# if __name__ == '__main__':
#     main()

