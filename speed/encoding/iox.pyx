# cdef char* alphabet = 'ACTG'
cdef char* alphabet = "ACTG"

def int2seq(int x, int base, int seqlen):
    cdef char *result =  malloc(seqlen*sizeof(char));
    cdef int modulo
    cdef int counter = 0
    if x == 0:
        return [alphabet[0]]


    while x:
        modulo  = x % base
        # print(alphabet[modulo])
        result[counter] = alphabet[modulo]
        x = x // base
    return result

#def itoa_wrapper(x, base):
#    target = ''
#    return itoa(x, target, 4)
