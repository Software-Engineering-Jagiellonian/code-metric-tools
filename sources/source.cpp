#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cstdarg>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

typedef std::pair<std::size_t, std::size_t> coordiates;
typedef std::vector<double> vector;
typedef std::vector<vector> matrix;
typedef std::function<coordiates(matrix, vector, std::size_t)> function;

coordiates maxPartialChoice(const matrix &A, const vector &b, std::size_t k){
    int n = b.size();
    coordiates max = coordiates(k,k);
    for(int i = k + 1 ; i < n ; i++){
        if(std::abs(A[i][k]) > std::abs(A[max.first][k]))
            max.first = i;
    }
    return max;
}

coordiates maxFullChoice(const matrix &A, const vector &b, std::size_t k){
    int n = b.size();
    coordiates max = coordiates(k, k);
    for(int i = k ; i < n ; i++){
        for(int j = k ; j < n ; j++){
            if(std::abs(A[i][j]) > std::abs(A[max.first][max.second]))
                max = coordiates(i, j);
        }
    }
    return max;
}

size_t rowMaxNorm(const vector &vec){
    int n = vec.size();
    double max = 0;
    for(int i = 1 ; i < n ; i++){
        if(std::abs(vec[i]) > std::abs(vec[max]))
            max = i;
    }
    return max;
}

coordiates maxScaledPartialChoice(const matrix &A, const vector &b, std::size_t k){
    int n = b.size();
    coordiates max = coordiates(k,k);
    for(int i = k + 1 ; i < n ; i++){
        if(std::abs(A[i][k] / A[i][rowMaxNorm(A[i])]) > std::abs(A[max.first][k] / A[max.first][rowMaxNorm(A[max.first])]))
            max.first = i;
    }
    return max;
}

void swap(matrix &mat, vector &vec, coordiates coord, int k){
    std::swap(mat[coord.first], mat[k]);
    if(coord.second != k){
        for(int i = 0 ; i < vec.size() ; i++){
            std::swap(mat[i][coord.second], mat[i][k]);
        }
        std::swap(vec[coord.second], vec[k]);
    }
}

void triangulation(matrix &mat, vector &vec, function const &method){
    int n = vec.size();
    for(int i = 0 ; i < n ; i++){
        coordiates maxCoord= method(mat, vec, i);
        swap(mat, vec, maxCoord, i); //zamienia odpowiednie kolumny i wierwsze, zwraca= coeff
        for(int j = i + 1 ; j < n ; j++){
            for(int k = i + 1 ; k < n ; k++){
                mat[j][k] -= (mat[j][i] / mat[i][i]) * mat[i][k];
            }
            vec[j] -= (mat[j][i] / mat[i][i]) * vec[i];
            mat[j][i] = 0;
        }
    }
}

vector Gauss(const matrix &A, const vector &b ,function const &method){
    std::size_t n = b.size();
    matrix mat = A;
    vector solution = b;
    triangulation(mat, solution, method);
    for(int i = n - 1 ; i >= 0 ; i--){
        for(int j = n - 1 ; j > i ; j--){
            solution[i] -= mat[i][j] * solution[j];
        }
        solution[i] = solution[i] / mat[i][i];
    }
    return solution;
}

void printVec(const vector &vec){
    for(auto i : vec)
        std::cout << i << " ";
}

void printMat(const matrix &mat){
    for(auto i : mat) {
        printVec(i);
        std::cout << std::endl;
    }
}

void readElem(char* str, int* elem, int count = 0){
    if(count == 5)
        return;
    else if(*str == '0'){
        (*elem) <<= 1;
    } else if(*str == '1') {
        (*elem) <<= 1;
        (*elem) += 1;
    }
    readElem(str + 1, elem, count + 1);
}

void readSet(char* str, int* set){
    if(*str == '\0')
        return;
    else if((*str == '0') || (*str == '1')){
        int elem = 0;
        readElem(str, &elem);
        int mask = (1 << elem);
        *set = (*set) | mask;
        str = str + 4;
    }
    readSet(str+1, set); ///<--- +1?

}
void Emplace(char* str, int* set){
    *set = 0;
    readSet(str, set);
}

void Insert(char* str, int* set){
    readSet(str, set);
}

void Erase(char* str, int* set){
    int tmp = 0;
    readSet(str, &tmp);
    *set = (*set) & (~tmp);
}

void mystrcpy(char* strl, const char* strr){
    *strl = *strr;
    if(*strr == '\0') return;
    mystrcpy(strl  + 1, strr + 1);
}

void printInElem(int elem, char* str, int mask = 16){
    if(mask == 0){
        return;
    }
    if((elem & mask) == mask){  //1
        *str = '1';
    } else { //0
        *str = '0';
    }
    printInElem(elem,str + 1, (mask >> 1));
}

void printInSet(int set, char* str, int first = 1, int elem = 31){
    if(elem == -1) {
        *str = '\0';
        return;
    }
    int mask = (1 << elem);
    if((set & mask) == mask){  //1 na bicie numer elem
        if(first == 1){
            printInElem(elem, str);
            first = 0;
            printInSet(set, str +  5, first, elem - 1);
        } else {
            *str = ' ';
            printInElem(elem, str + 1);
            printInSet(set, str +  6, first, elem - 1);
        }
    } else { //0
        printInSet(set, str, first, elem - 1);
    }
}

void Print(int set, char* str){
    if(set == 0){  // set = 0 <-> pusty zbior
        strcpy(str, "empty");
        return;
    }
    printInSet(set, str);

}

bool Emptiness(int set){
    return (set == 0);
}
bool Nonempty(int set){
    return !Emptiness(set);
}
bool Member(char* str, int set){
    int tmp = 0;
    readSet(str, &tmp);
    return ((set & tmp) == tmp);
}
bool Disjoint(int set1, int set2){
    return ((set1 & set2) == 0);
}
bool Conjunctive(int set1, int set2){
    return !Disjoint(set1, set2);
}
bool Equality(int set1, int set2){
    return (set1 == set2);
}
bool Inclusion(int set1, int set2){
    return ((set2 & set1) == set1);
}
void Union(int set1, int set2, int* unionset){
    *unionset = (set1 | set2);
}
void Intersection(int set1, int set2, int* intersection){
    *intersection = (set1 & set2);
}
void Symmetric(int set1, int set2, int* symdiff){
    int unionset;
    Union(set1, set2, &unionset);
    int intersection;
    Intersection(set1, set2, &intersection);
    *symdiff = unionset & (~intersection);
}
void Difference(int set1, int set2, int* diff){
    int intersection;
    Intersection(set1, set2, &intersection);
    *diff = set1 & (~intersection);
}
void Complement(int set, int* complement){
    *complement = (~set);
}

int sumBytes(int set, int mask = 1){
    if(mask == 0)
        return 0;
    else if((mask & set) == mask)
        return 1 + sumBytes(set, mask << 1);
    else
        return 0 + sumBytes(set, mask << 1);
}

int Cardinality(int set){
    return sumBytes(set);
}

bool LessThen(int set1, int set2){
    int set1card = Cardinality(set1);
    int set2card = Cardinality(set2);
    if(set1card != set2card)
        return (set1card < set2card);
    return set1 < set2;
}
bool GreatThen(int set1, int set2){
    int set1card = Cardinality(set1);
    int set2card = Cardinality(set2);
    if(set1card != set2card)
        return (set1card > set2card);
    return set1 > set2;
}

bool LessEqual(int set1, int set2){
    return !GreatThen(set1, set2);
}
bool GreatEqual(int set1, int set2){
    return !LessThen(set1, set2);
}

struct Node{
    int id;
    char id1, id2, id3;
    Node* next;
    Node* prev;
    Node* subList;

    Node(char id1, char id2, char id3, int id = -1){
        this->id = id;
        this->id1 = id1;
        this->id2 = id2;
        this->id3 = id3;
        next = NULL;
        prev = NULL;
        subList = NULL;
    }
};
//--------------------------------------------------------------------------------------------------------------------//
//jezeli d1 > d2 zwraca 1, 0  gdy rowne, -1 gdy d1 < d2
int idcmp(char d11, char d12, char d13, char d21, char d22, char d23){
    if(d11 == d21 && d12 == d22 && d13 == d23) return 0;
    if(d11 > d21) return 1;
    if(d11 == d21 && d12 > d22) return 1;
    if(d11 == d21 && d12 == d22 && d13 > d23) return 1;
    return -1;
}
//jezeli lNode > rNode zwraca 1, jak rowne zwraca 0, -1 gdzy lNode < rNode
int nodecmp(Node* lNode, Node* rNode){
    int cmp = idcmp(lNode->id1, lNode->id2, lNode->id3, rNode->id1, rNode->id2, rNode->id3);
    if(cmp > 0)
        return 1;
    else if(cmp == 0){
        if(lNode->id > rNode->id)
            return 1;
        else if(lNode->id < rNode->id)
            return -1;
        else
            return 0;
    }
    return -1;
}
//--------------------------------------------------------------------------------------------------------------------//
//dodanie elementu do listy, elementy dodawane sa z zachowaniem porzadku z zadania
void addNode(Node** list, Node* newElem){
    //dodanie do pustej listy
    if(*list == NULL){
        *list = newElem;
        return;
    }
    Node* tmp = *list;
    while(tmp->next != NULL) {
        int compare = nodecmp(newElem, tmp);
        if (compare < 0) {
            if (tmp->prev == NULL)
                *list = newElem;
            else {
                tmp->prev->next = newElem;
                newElem->prev = tmp->prev;
            }
            newElem->next = tmp;
            tmp->prev = newElem;
            return;
        }
        tmp = tmp->next;
    }
    int compare = nodecmp(newElem, tmp);
    if((tmp->prev == NULL) && (compare < 0)){
        tmp->prev = newElem;
        newElem->next = tmp;
        *list =newElem;
        return;
    }
    if(compare < 0){
        tmp->prev->next = newElem;
        newElem->prev = tmp->prev;
        newElem->next = tmp;
        tmp->prev = newElem;
        return;
    }
    //dodajemy na koncu
    newElem->prev = tmp;
    tmp->next = newElem;
}
//--------------------------------------------------------------------------------------------------------------------//
//usuwanie zaalokowanej dynamicznie pamieci, funkcja wywoluje sie rekurencyjnie na wszystkich podlistach
void clearList(Node** list){
    Node* delElem = *list;
    while(delElem != NULL){
        Node* tmp = delElem->next;
        clearList(&(delElem->subList));
        delete delElem;
        delElem = tmp;
    }
    *list = NULL;
}
//--------------------------------------------------------------------------------------------------------------------//
//funkcja zwraca adres wskaznika na podliste w zadanym wezle listy
Node** findSublist(Node** list, Node node){
    Node* tmp = *list;
    while(tmp != NULL){
        if(nodecmp(tmp, &node) == 0)
            return &(tmp->subList);
        tmp = tmp->next;
    }
    return NULL;
}
//--------------------------------------------------------------------------------------------------------------------//
//zwraca wskaznik na szukany wezel, jak nie znajdzie zwraca NULL
Node* findNode(Node** list, Node node){
    Node* tmp = *list;
    while(tmp != NULL){
        if(nodecmp(tmp, &node) == 0)
            return tmp;
        tmp = tmp->next;
    }
    return NULL;
}
//--------------------------------------------------------------------------------------------------------------------//
//usuwa wezel
void deleteNode(Node** list, Node* node){
    if(node == NULL) return;
    if(*list == NULL) return;
    if(node->prev == NULL){
        *list = node->next;
        if(*list != NULL)
            node->next->prev = NULL;
    } else if(node->next == NULL) {
        node->prev->next = NULL;
    } else {
        node->prev->next = node->next;
        node->next->prev = node->prev;
    }
    delete node;
}
//--------------------------------------------------------------------------------------------------------------------//
//operacja P, trzy funkcje do wyswietlania list
void printSubsublist(Node** subsublist){
    Node* tmp = *subsublist;
    while(tmp != NULL){
        std::cout << " " << tmp->id1 << tmp->id2 << tmp->id3 << " " << tmp->id;
        tmp = tmp->next;
    }
}
void printSublist(Node** sublist){
    Node* tmp = *sublist;
    while(tmp != NULL){
        std::cout << tmp->id;
        printSubsublist(&(tmp->subList));
        std::cout << std::endl;
        tmp = tmp->next;
    }
}
void printList(Node** list){
    Node* tmp = *list;
    while(tmp != NULL){
        std::cout << tmp->id1 << tmp->id2 << tmp->id3 << std::endl;
        printSublist(&(tmp->subList));
        tmp = tmp->next;
    }
}
//--------------------------------------------------------------------------------------------------------------------//
//funkcja pomocnicza do tworzenia nowego wezla listy
Node* newNode(char c1, char c2, char c3, int id = -1){
    Node* newElem = new Node(c1, c2, c3, id);
    return newElem;
}
//--------------------------------------------------------------------------------------------------------------------//
//Usuwa wszystie wezly w liscie o takich identyfikatorac jak w node
void deleteMatchingNodes(Node ** list, Node node){
    Node* listtmp = *list;
    while(listtmp != NULL){
        Node* sublisttmp = listtmp->subList;
        while (sublisttmp != NULL){
            Node* delElem = findNode(&(sublisttmp->subList), node);
            if((delElem != NULL) && (nodecmp(delElem, &node) == 0))
                deleteNode(&(sublisttmp->subList), delElem);
            sublisttmp = sublisttmp->next;
        }
        listtmp = listtmp->next;
    }
}
//--------------------------------------------------------------------------------------------------------------------//
//przenosi podany wezel na podana liste
void moveNode(Node** from, Node** to, Node* node){
    if(node == NULL) return;
    if(node->prev == NULL){
        *from = node->next;
        if(*from != NULL)
            node->next->prev = NULL;
    } else if(node->next == NULL) {
        node->prev->next = NULL;
    } else {
        node->prev->next = node->next;
        node->next->prev = node->prev;
    }
    node->prev = NULL;
    node->next = NULL;
    addNode(to, node);
}
//--------------------------------------------------------------------------------------------------------------------//
//operacja R
void deleteNodeRecursive(Node** list, Node* node){
    Node* delNode = node->subList;
    node->subList = NULL;
    while(delNode != NULL){
        clearList(&(delNode->subList));
        deleteMatchingNodes(list, *delNode);
        Node* tmp = delNode->next;
        delete delNode;
        delNode = tmp;
    }
    deleteNode(list, node);
}
//--------------------------------------------------------------------------------------------------------------------//
//funkcja wywolujaca komendy
void callCommand(char command, Node** list){
    switch (command) {
        case 'S': {
            char c1, c2, c3;
            std::cin >> c1 >> c2 >> c3;
            addNode(list, newNode(c1, c2, c3));
            break;
        }
        case 'P': {
            printList(list);
            break;
        }
        case 'B': {
            char c1, c2, c3;
            int id;
            std::cin >> id >> c1 >> c2 >> c3;
            addNode(findSublist(list, Node(c1, c2, c3)), newNode(c1, c2, c3, id));
            break;
        }
        case 'L': {
            char c11, c12, c13, c21, c22, c23;
            int id1, id2;
            std::cin >> id1 >> c11 >> c12 >> c13 >> id2 >> c21 >> c22 >> c23;
            addNode(findSublist(findSublist(list, Node(c11, c12, c13)),Node(c11, c12, c13, id1)),
                    newNode(c21, c22, c23, id2));
            break;
        }
        case 'U': {
            char c11, c12, c13, c21, c22, c23;
            int id1, id2;
            std::cin >> id1 >> c11 >> c12 >> c13 >> id2 >> c21 >> c22 >> c23;
            Node** subsublist = findSublist(findSublist(list, Node(c11, c12, c13)),
                                            Node(c11, c12, c13, id1));
            deleteNode(subsublist,findNode(subsublist, Node(c21, c22, c23, id2)));
            break;
        }
        case 'D': {
            char c1, c2, c3;
            int id;
            std::cin >> id >> c1 >> c2 >> c3;
            Node** sublist = findSublist(list, Node(c1, c2, c3));
            Node** subsublist = findSublist(sublist, Node(c1, c2, c3, id));
            clearList(subsublist);
            deleteNode(sublist, findNode(sublist, Node(c1, c2, c3, id)));
            deleteMatchingNodes(list, Node(c1, c2, c3, id));
            break;
        }
        case 'M': {
            char c11, c12, c13, c21, c22, c23;
            int id;
            std::cin >> id >> c11 >> c12 >> c13 >> c21 >> c22 >> c23;
            Node** from = findSublist(list, Node(c11, c12, c13));
            Node* node = findNode(from,Node(c11, c12, c13, id));
            node->id1 = c21;
            node->id2 = c22;
            node->id3 = c23;
            Node** to = findSublist(list, Node(c21, c22, c23));
            moveNode(from, to, node);
            break;
        }
        case 'R': {
            char c1, c2, c3;
            std::cin >> c1 >> c2 >> c3;
            deleteNodeRecursive(list, findNode(list, Node(c1, c2, c3)));
            break;
        }
    }
}
//--------------------------------------------------------------------------------------------------------------------//
void makeArrayFromArgs(std::string* argsArray, int numOfArgs, va_list args){
    if(numOfArgs == 0)
        return;
    argsArray[numOfArgs - 1] = std::string(va_arg(args, char*));
    makeArrayFromArgs(argsArray,numOfArgs - 1, args);
}


//return numerical vaulue of char (if char is not in [0-9] returns 0)
int numericalVal(char sign){
    int val = sign - '0';
    if(val >= 0 && val <= 9) return val;
    return 0;
}

//std::string cutFront(const std::string &number, int idx){
std::string cutFront(const std::string number, int idx = 0){
    if(numericalVal(number[idx]) > 0){
        return number.substr(idx);
    }
    if(idx == number.length()) return "0";
    return cutFront(number, idx + 1);
}

//sign of number in string
std::string sign(const std::string &number){
    if(number[0] == '-') return "-";
    return "";
}

//substracts two numbers wrtten in strings, returns string with difference, ignores "+" or "-" at the beginning
std::string substractRecursive(const std::string& greatertArg, const std::string& smallerArg, int carryValue = 0, int idxFromLast = 1){
//std::string substractRecursive(const std::string& greatertArg, const std::string& smallerArg, int carryValue, int idxFromLast){
    if(idxFromLast > greatertArg.length()){
        return "";
    }
    int sum;
    if(idxFromLast <= smallerArg.length()){
        sum = numericalVal(greatertArg[greatertArg.length() - idxFromLast]) - numericalVal(smallerArg[smallerArg.length() - idxFromLast]) + carryValue;
    } else{ // if(idxFromLast <= greatertArg.length()){
        sum = numericalVal(greatertArg[greatertArg.length() - idxFromLast])  + carryValue;
    }
    carryValue = 0;
    if(sum < 0){
        sum += 10;
        carryValue = -1;
    }
    return (substractRecursive(greatertArg, smallerArg, carryValue, idxFromLast + 1) + (char)(sum + '0'));
}


//adds recursive two numbers written in strings, returns string with sum, ignores "+" or "-" at the beginning
std::string addRecursive(const std::string& firstArg, const std::string& secondArg, int carryValue = 0, int idxFromLast = 1){
//std::string addRecursive(const std::string& firstArg, const std::string& secondArg, int carryValue, int idxFromLast){
    if(idxFromLast > firstArg.length() && idxFromLast > secondArg.length()){
        if(carryValue > 0){
            std::string str;
            return (char)(carryValue + '0') + str;
        }
        return "";
    }
    int sum;
    if(idxFromLast <= firstArg.length() && idxFromLast <= secondArg.length()){
        sum = numericalVal(firstArg[firstArg.length() - idxFromLast]) + numericalVal(secondArg[secondArg.length() - idxFromLast]) + carryValue;
    } else if(idxFromLast <= firstArg.length()){
        sum = numericalVal(firstArg[firstArg.length() - idxFromLast]) + carryValue;
    } else{ //if(idxFromLast <= secondArg.length()){
        sum = numericalVal(secondArg[secondArg.length() - idxFromLast]) + carryValue;
    }
    carryValue = sum / 10;
    return (addRecursive(firstArg, secondArg, carryValue, idxFromLast + 1) + (char)((sum % 10) + '0'));
}

//returns sum of n numbers written in strings
std::string Sum(int n, const std::string* numbers){
    if(n == 2){
        std::string arg1 = cutFront(numbers[0]);
        std::string arg2 = cutFront(numbers[1]);
        std::string sum;
        if(sign(numbers[0]) == sign(numbers[1])){
            sum = sign(numbers[0]) + cutFront(addRecursive(arg1, arg2));
        } else{
            if(arg1.compare(arg2) > 0)  //|numbers[0]| > |numbers[1]|
                sum = sign(numbers[0]) + cutFront(substractRecursive(arg1, arg2));
            else
                sum = sign(numbers[1]) + cutFront(substractRecursive(arg2, arg1));
        }
        if(cutFront(sum) == "0") return "0";
        return sum;
    } else{
        std::string prevSum = Sum(n - 1, numbers);
        std::string arg1 = cutFront(prevSum);
        std::string arg2 = cutFront(numbers[n - 1]);
        std::string sum;
        if(sign(prevSum) == sign(numbers[n - 1])){
            sum = sign(prevSum) + cutFront(addRecursive(arg1, arg2));
        } else{
            if(arg1.compare(arg2) > 0)
                sum = sign(prevSum) + cutFront(substractRecursive(arg1, arg2));
            else
                sum = sign(numbers[n - 1]) + cutFront(substractRecursive(arg2, arg1));
        }
        if(cutFront(sum) == "0") return "0";
        return sum;
    }
}

void Sum(std::string* sum, int n, const std::string* numbers){
    *sum = Sum(n, numbers);
}

void Sum(std::string& sum, int n, const std::string* numbers){
    sum = Sum(n, numbers);
}

std::string vSum(int n, va_list vNumbers){
    if(n == 2){
        std::string arg1 = (std::string)(va_arg(vNumbers, char*));
        std::string arg2 = (std::string)(va_arg(vNumbers, char*));
        std::string sum;
        if(sign(arg1) == sign(arg2)){
            sum = sign(arg1) + cutFront(addRecursive(arg1, arg2));
        } else{
            if(cutFront(arg1).compare(cutFront(arg2)) > 0)
                sum = sign(arg1) + cutFront(substractRecursive(arg1, arg2));
            else
                sum = sign(arg2) + cutFront(substractRecursive(arg2, arg1));
        }
        if(cutFront(sum) == "0") return "0";
        return sum;
    } else{
        std::string arg1 = vSum(n - 1, vNumbers);
        std::string arg2 = (std::string)(va_arg(vNumbers, char*));
        std::string sum;
        if(sign(arg1) == sign(arg2)){
            sum = sign(arg1) + cutFront(addRecursive(arg1, arg2));
        } else{
            if(cutFront(arg1).compare(cutFront(arg2)) > 0)
                sum = sign(arg1) + cutFront(substractRecursive(arg1, arg2));
            else
                sum = sign(arg2) + cutFront(substractRecursive(arg2, arg1));
        }
        if(cutFront(sum) == "0") return "0";
        return sum;
    }
}

std::string Sum(int n, ...){
    std::string sum;
    va_list args;
    va_start(args, n);
    sum = vSum(n, args);
    va_end(args);
    return sum;
}

void Sum(std::string* sum, int n, ...){
    std::string str;
    std::string argsArr[n];
    va_list args;
    va_start(args, n);
    makeArrayFromArgs(argsArr, n, args);
    //str = vSum(n, args);
    str = Sum(n, argsArr);
    va_end(args);
    (*sum).erase();
    *sum = str;
}

void Sum(std::string& sum, int n, ...){
    std::string str;
    std::string argsArr[n];
    sum.erase();
    va_list args;
    va_start(args, n);
    makeArrayFromArgs(argsArr, n, args);
    //sum = vSum(n, args);
    sum = Sum(n, argsArr);
    va_end(args);
}

//-----------------------------------------------------------------------------------------------------------------------


//std::string multiplyByConstantRec(const std::string &number, const char constant, int idxFromLast, int carry){
std::string multiplyByConstantRec(const std::string &number, char &constant, int idxFromLast = 1, int carry = 0){
    if(idxFromLast > number.length()){
        std::string str;
        if(carry > 0)
            return (char)(carry + '0') + str;
        return str;
    }
    int a = numericalVal(number[number.length() - idxFromLast]) * numericalVal(constant) + carry;
    return multiplyByConstantRec(number, constant, idxFromLast + 1, a / 10) + (char)((a % 10) + '0');
}

//std::string multiplyRecursive(const std::string &firstArg, const std::string &secondArg, int idx){
std::string multiplyRecursive(std::string firstArg, std::string secondArg, int idx = 1){
    if((firstArg == "0") || (secondArg== "0")) return "0";
    if(firstArg == "1") return secondArg;
    if(secondArg == "1") return firstArg;
    if((firstArg.length() - idx) == 0) {
        return multiplyByConstantRec(secondArg, firstArg[0]);
    }
    std::string firstSumArg = multiplyByConstantRec(secondArg, firstArg[firstArg.length() - idx]);
    std::string secondSumArg = multiplyRecursive(firstArg, secondArg, idx + 1) + "0";
    return addRecursive(firstSumArg,secondSumArg);
}

std::string Mult(int n, const std::string* numbers){
    if(n == 2){
        std::string str = multiplyRecursive(cutFront(numbers[0]), cutFront(numbers[1]));
        if(cutFront(str) == "0") return "0";
        if(sign(numbers[0]) != sign(numbers[1]))
            return "-" + str;
        return str;
    } else{
        std::string prevMult = Mult(n - 1, numbers);
        std::string mult = multiplyRecursive(cutFront(prevMult), cutFront(numbers[n - 1]));
        if(mult == "0") return "0";
        if((sign(prevMult)) != (sign(numbers[n - 1]))){
            return "-" + mult;
        }
        return mult;
    }
}

std::string vMult(int n, va_list vNumbers){
    if(n == 2){
        std::string arg1 = (std::string)(va_arg(vNumbers, char*));
        std::string arg2 = (std::string)(va_arg(vNumbers, char*));
        std::string str = multiplyRecursive(cutFront(arg1), cutFront(arg2));
        if(cutFront(str) == "0") return "0";
        if(sign(arg1) != sign(arg2))
            return "-" + str;
        return str;
    } else{
        std::string prevMult = vMult(n - 1, vNumbers);
        std::string arg = (std::string)(va_arg(vNumbers, char*));
        std::string mult = multiplyRecursive(cutFront(prevMult), cutFront(arg));
        if(mult == "0") return "0";
        if((sign(prevMult)) != (sign(arg))){
            return "-" + mult;
        }
        return mult;
    }
}

void Mult(std::string* mult, int n, const std::string* numbers){
    *mult = Mult(n, numbers);
}

void Mult(std::string &mult, int n, const std::string* numbers){
    mult = Mult(n, numbers);
}


std::string Mult(int n, ...){
    std::string mult;
    va_list args;
    va_start(args, n);
    mult = vMult(n, args);
    va_end(args);
    return mult;
}

void Mult(std::string* mult, int n, ...){
    std::string argsArr[n];
    va_list args;
    va_start(args, n);
    makeArrayFromArgs(argsArr, n, args);
    //*mult = vMult(n, argsArr);
    *mult = Mult(n, argsArr);
    va_end(args);
}

void Mult(std::string& mult, int n, ...){
    std::string argsArr[n];
    va_list args;
    va_start(args, n);
    makeArrayFromArgs(argsArr, n, args);
    //mult = vMult(n, argsArr);
    mult = Mult(n, argsArr);
    va_end(args);
}
//-----------------------------------------------------------------------------------------------------------------------

std::string Operation(std::string(*operation)(int, const std::string*), int n, const std::string* numbers){
    return operation(n, numbers);
}

void Operation(std::string* result, std::string(* operation)(int, const std::string* ), int n, const std::string* numbers){
    *result = operation(n, numbers);

}

void Operation(std::string& result, void(* operation)(std::string*, int, const std::string* ), int n, const std::string* numbers){
    operation(&result, n, numbers);
}

std::string Operation(std::string(* operation)(int, const std::string*), int n, ...){
    std::string str;
    std::string argsArr[n];
    va_list args;
    va_start(args, n);
    makeArrayFromArgs(argsArr, n, args);
    str = operation(n, argsArr);
    va_end(args);
    return str;
}

void Operation(std::string* result, std::string(* operation)(int, const std::string* ), int n, ...){
    std::string argsArr[n];
    va_list args;
    va_start(args, n);
    makeArrayFromArgs(argsArr, n, args);
    *result = operation(n, argsArr);
    va_end(args);
}

void Operation(std::string& result, void(* operation)(std::string*, int, const std::string* ), int n, ...){
    std::string argsArr[n];
    va_list args;
    va_start(args, n);
    makeArrayFromArgs(argsArr, n, args);
    operation(&result, n, argsArr);
    va_end(args);
}

class FRUIT_CLASS;
class BRANCH_CLASS;
class TREE_CLASS;
class GARDEN_CLASS;

class FRUIT_CLASS{
private:
    unsigned int weight;
    unsigned int length;
    BRANCH_CLASS *branch;
    FRUIT_CLASS *nextFruit;

public:
    FRUIT_CLASS() : weight(0), length(0), branch(NULL), nextFruit(NULL) {};
    FRUIT_CLASS(unsigned int length, BRANCH_CLASS* branch)
            : weight(0), length(length), branch(branch), nextFruit(NULL) {};
    FRUIT_CLASS(const FRUIT_CLASS &obj)
            : weight(obj.getWeight()), length(obj.getLength()), branch(NULL), nextFruit(NULL) {};

    FRUIT_CLASS* getNext() { return nextFruit; }
    void setNext(FRUIT_CLASS* fruit) { this->nextFruit = fruit; }

    unsigned int getLength() const { return length; }
    unsigned int getWeight() const { return weight; }
    void growthFruit();
    void fadeFruit();
    void pluckFruit();
    BRANCH_CLASS* getBranchPointer() const { return branch; }

    void setBranch(BRANCH_CLASS *obj) { branch = obj; }

    void print() const;
};

class BRANCH_CLASS{
private:
    unsigned int fruitsCount;
    unsigned int fruitsWeight;
    unsigned int height;
    unsigned int length;
    BRANCH_CLASS *nextBranch;
    FRUIT_CLASS* fruits;
    TREE_CLASS* tree;
public:
    BRANCH_CLASS()
            : fruitsCount(0), fruitsWeight(0), height(0), length(0), nextBranch(NULL),  fruits(NULL), tree(NULL) {};
    BRANCH_CLASS(unsigned int height, TREE_CLASS *tree)
            : fruitsCount(0), fruitsWeight(0), height(height), length(0), nextBranch(NULL),  fruits(NULL),
              tree(tree) {};
    BRANCH_CLASS(const BRANCH_CLASS &);
    ~BRANCH_CLASS();

    BRANCH_CLASS* getNext() { return nextBranch; }
    void setNext(BRANCH_CLASS *branch) { nextBranch = branch; }

    unsigned int getFruitsTotal() const { return fruitsCount; }
    unsigned int getWeightsTotal() const { return fruitsWeight; }
    unsigned int getHeight() const { return height; }
    unsigned int getLength() const { return length; }
    void growthBranch();
    void fadeBranch();
    void harvestBranch(unsigned int);
    void cutBranch(unsigned int);
    FRUIT_CLASS* getFruitPointer(unsigned int);
    TREE_CLASS* getTreePointer() const { return tree; }

    FRUIT_CLASS* getFruitList() const { return fruits; }
    void setTree(TREE_CLASS *tree) { this->tree = tree; }
    void setHeight(unsigned int height) { this->height = height; }
    void setWeight(unsigned int weight) { fruitsWeight = weight; }

    void print() const;
};

class TREE_CLASS{
private:
    unsigned int id;
    unsigned int height;
    unsigned int branchesCount;
    unsigned int fruitsCount;
    unsigned int fruitsWeight;
    BRANCH_CLASS *branches;
    GARDEN_CLASS *garden;
    TREE_CLASS *nextTree;
    TREE_CLASS *prevTree;
public:
    TREE_CLASS()
            : id(0), height(0), branchesCount(0), fruitsCount(0), fruitsWeight(0), branches(NULL), garden(NULL),
              nextTree(NULL), prevTree(NULL) {};
    TREE_CLASS(unsigned int id, GARDEN_CLASS *garden)
            : id(id), height(0), branchesCount(0), fruitsCount(0), fruitsWeight(0), branches(NULL), garden(garden),
              nextTree(NULL), prevTree(NULL) {};
    TREE_CLASS(const TREE_CLASS &);
    ~TREE_CLASS();

    void setNext(TREE_CLASS *tree) { nextTree = tree; }
    TREE_CLASS* getNext() { return nextTree; }
    void setPrev(TREE_CLASS *tree) { prevTree = tree; }
    TREE_CLASS* getPrev() { return prevTree; }

    unsigned int getBranchesTotal() const { return branchesCount; }
    unsigned int getFruitsTotal() const { return fruitsCount; }
    unsigned int getWeightsTotal() const { return fruitsWeight; }
    unsigned int getNumber() const { return id; }
    unsigned int getHeight() const { return height; }
    void growthTree();
    void fadeTree();
    void harvestTree(unsigned int);
    void cutTree(unsigned int);
    void cloneBranch(BRANCH_CLASS*);
    GARDEN_CLASS* getGardenPointer() const { return garden; }
    BRANCH_CLASS* getBranchPointer(unsigned int);

    BRANCH_CLASS* getBranchesList() const { return branches; }
    void setNumber(unsigned int num) { id = num; }
    void setGarden(GARDEN_CLASS *garden) { this->garden = garden; }
    void setWeight(unsigned int weight) { fruitsWeight = weight; }
    void setFruitCount(unsigned int count) { fruitsCount = count; }


    void print();
};

class GARDEN_CLASS{
private:
    unsigned int treesCount;
    unsigned int branchesCount;
    unsigned int fruitsCount;
    unsigned int fruitsWeight;
    TREE_CLASS *trees;
    TREE_CLASS *lastTree;
public:
    GARDEN_CLASS() : treesCount(0), branchesCount(0), fruitsCount(0), fruitsWeight(0), trees(NULL), lastTree(NULL) {};
    ~GARDEN_CLASS();

    unsigned int getTreesTotal() { return treesCount; }
    unsigned int getBranchesTotal() { return branchesCount; }
    unsigned int getFruitsTotal() { return fruitsCount; }
    unsigned int getWeightsTotal() { return fruitsWeight; }
    void plantTree();
    void extractTree(unsigned int);
    void growthGarden();
    void fadeGarden();
    void harvestGarden(unsigned int);
    TREE_CLASS* getTreePointer(unsigned int);
    void cloneTree(unsigned int);

    void setWeight(unsigned int weight) { fruitsWeight = weight; }
    void setFruitCount(unsigned int count) { fruitsCount = count; }
    void setBranchesCount(unsigned int count) { branchesCount = count; }

    void print();
};


void FRUIT_CLASS::print() const{
    std::cout << "   > fruit at length " << length << ": weight = " << weight << std::endl;
}

void BRANCH_CLASS::print() const{
    std::cout << "  >> branch at height " << height << ": length = " << length << ", fruits count = " << fruitsCount
              << ", fruits weight = " << fruitsWeight << std::endl;
    FRUIT_CLASS *fruit = fruits;
    while(fruit != NULL){
        fruit->print();
        fruit = fruit->getNext();
    }
}

void TREE_CLASS::print() {
    std::cout << " >>> tree with id " << id << ": height = " << height << ", branches count = " << branchesCount
              << " fruits count = " << fruitsCount << ", fruits weight = " << fruitsWeight << std::endl;
    BRANCH_CLASS *branch = branches;
    while(branch != NULL){
        branch->print();
        branch = branch->getNext();
    }
}

void GARDEN_CLASS::print() {
    std::cout << ">>>> garden: trees count = " << treesCount << ", branches count = " << branchesCount
              << ", fruits count = " << fruitsCount << ", fruits weight = " << fruitsWeight << std::endl;
    TREE_CLASS *tree = trees;
    while(tree != NULL){
        tree->print();
        tree = tree->getNext();
    }
}

//--- FRUIT_CLASS ------------------------------------------------------------------------------------------------------

void FRUIT_CLASS::growthFruit() {
    weight++;
    if(branch != NULL){
        branch->setWeight(branch->getWeightsTotal() + 1);
        TREE_CLASS *tree = branch->getTreePointer();
        if(tree != NULL){
            tree->setWeight(tree->getWeightsTotal() + 1);
            GARDEN_CLASS *garden = tree->getGardenPointer();
            if(garden != NULL)
                garden->setWeight(garden->getWeightsTotal() + 1);
        }
    }
}

void FRUIT_CLASS::fadeFruit() {
    if(weight == 0) return;
    weight--;
    if(branch != NULL){
        branch->setWeight(branch->getWeightsTotal() - 1);
        TREE_CLASS *tree = branch->getTreePointer();
        if(tree != NULL){
            tree->setWeight(tree->getWeightsTotal() - 1);
            GARDEN_CLASS *garden = tree->getGardenPointer();
            if(garden != NULL)
                garden->setWeight(garden->getWeightsTotal() - 1);
        }
    }
}

void FRUIT_CLASS::pluckFruit() {
    if(weight == 0) return;
    if(branch != NULL){
        branch->setWeight(branch->getWeightsTotal() - weight);
        TREE_CLASS *tree = branch->getTreePointer();
        if(tree != NULL){
            tree->setWeight(tree->getWeightsTotal() - weight);
            GARDEN_CLASS *garden = tree->getGardenPointer();
            if(garden != NULL)
                garden->setWeight(garden->getWeightsTotal() - weight);
        }
    }
    weight = 0;
}

//--- BRANCH_CLASS -----------------------------------------------------------------------------------------------------

BRANCH_CLASS::BRANCH_CLASS(const BRANCH_CLASS &branch) {
    height = branch.getHeight();
    length = branch.getLength();
    fruitsCount = branch.getFruitsTotal();
    fruitsWeight = branch.getWeightsTotal();
    tree = NULL;
    nextBranch = NULL;

    fruits = NULL;
    FRUIT_CLASS *fruitToCopy = branch.getFruitList();

    if(fruitToCopy != NULL) {
        fruits = new FRUIT_CLASS(*fruitToCopy);
        fruits->setBranch(this);
        fruitToCopy = fruitToCopy->getNext();

        FRUIT_CLASS *fruit = fruits;
        while (fruitToCopy != NULL) {
            FRUIT_CLASS *clone = new FRUIT_CLASS(*fruitToCopy);
            clone->setBranch(this);
            fruit->setNext(clone);
            fruit = fruit->getNext();
            fruitToCopy = fruitToCopy->getNext();
        }
    }
}

BRANCH_CLASS::~BRANCH_CLASS() {
    while(fruits != NULL){
        FRUIT_CLASS *toDelete = fruits;
        fruits = fruits->getNext();
        delete toDelete;
    }
}

void BRANCH_CLASS::growthBranch() {
    length++;
    FRUIT_CLASS *fruit = fruits;
    while(fruit != NULL){
        fruit->growthFruit();
        fruit = fruit->getNext();
    }
    if(length % 2 == 0){
        FRUIT_CLASS *newFruit = new FRUIT_CLASS(length, this);
        if(fruitsCount > 0)
            newFruit->setNext(fruits);
        fruits = newFruit;
        fruitsCount++;
        if(tree != NULL){
            tree->setFruitCount(tree->getFruitsTotal() + 1);
            GARDEN_CLASS *garden = tree->getGardenPointer();
            if(garden != NULL)
                garden->setFruitCount(garden->getFruitsTotal() + 1);
        }
    }
}

void BRANCH_CLASS::fadeBranch() {
    if(length == 0) return;
    FRUIT_CLASS *fruit = fruits;
    while(fruit != NULL){
        fruit->fadeFruit();
        fruit = fruit->getNext();
    }
    if(fruits != NULL){
        if(fruits->getLength() == length){
            FRUIT_CLASS *toDelete = fruits;
            toDelete->pluckFruit();
            fruits = fruits->getNext();
            fruitsCount--;
            delete toDelete;
            if(tree != NULL){
                tree->setFruitCount(tree->getFruitsTotal() - 1);
                GARDEN_CLASS *garden = tree->getGardenPointer();
                if(garden != NULL)
                    garden->setFruitCount(garden->getFruitsTotal() - 1);
            }
        }
    }
    length--;
}

void BRANCH_CLASS::harvestBranch(unsigned int weight) {
    FRUIT_CLASS *fruit = fruits;
    while(fruit != NULL){
        if(fruit->getWeight() >= weight) {
            fruit->pluckFruit();
        }
        fruit = fruit->getNext();
    }
}

void BRANCH_CLASS::cutBranch(unsigned int length) {
    if(length >= this->length) return;
    this->length = length;
    if(fruits == NULL) return;
    FRUIT_CLASS *fruit = fruits;
    while(fruit != NULL){
        if(fruit->getLength() > length){
            FRUIT_CLASS *toDelete = fruit;
            toDelete->pluckFruit();
            fruits = fruits->getNext();
            fruitsCount--;
            if(tree != NULL){
                tree->setFruitCount(tree->getFruitsTotal() - 1);
                GARDEN_CLASS *garden = tree->getGardenPointer();
                if(garden != NULL)
                    garden->setFruitCount(garden->getFruitsTotal() - 1);
            }
            fruit = fruit->getNext();
            delete toDelete;
        } else
            fruit = fruit->getNext();
    }
}

FRUIT_CLASS* BRANCH_CLASS::getFruitPointer(unsigned int length) {
    FRUIT_CLASS *fruit = fruits;
    while(fruit != NULL){
        if(fruit->getLength() == length) return fruit;
        fruit = fruit->getNext();
    }
    return fruit;
}

//--- TREE_CLASS -------------------------------------------------------------------------------------------------------

TREE_CLASS::TREE_CLASS(const TREE_CLASS &obj) {
    id = 0;
    height = obj.getHeight();
    branchesCount = obj.getBranchesTotal();
    fruitsCount = obj.getFruitsTotal();
    fruitsWeight = obj.getWeightsTotal();
    garden = NULL;
    nextTree = NULL;
    prevTree = NULL;

    branches = NULL;
    BRANCH_CLASS *branchToCLone = obj.getBranchesList();

    if(branchToCLone != NULL){
        branches = new BRANCH_CLASS(*branchToCLone);
        branches->setTree(this);
        branchToCLone = branchToCLone->getNext();

        BRANCH_CLASS *branch = branches;
        while(branchToCLone != NULL){
            BRANCH_CLASS *clone = new BRANCH_CLASS(*branchToCLone);
            clone->setTree(this);
            branch->setNext(clone);
            branch = branch->getNext();
            branchToCLone = branchToCLone->getNext();
        }
    }
}

TREE_CLASS::~TREE_CLASS() {
    while(branches != NULL){
        BRANCH_CLASS *toDelete = branches;
        branches = branches->getNext();
        delete toDelete;
    }
}

void TREE_CLASS::growthTree() {
    height++;
    BRANCH_CLASS *branch = branches;
    while(branch != NULL){
        branch->growthBranch();
        branch = branch->getNext();
    }
    if(height % 3 == 0){
        BRANCH_CLASS *newBranch = new BRANCH_CLASS(height, this);
        if(branchesCount > 0)
            newBranch->setNext(branches);
        branches = newBranch;
        branchesCount++;
        if(garden != NULL)
            garden->setBranchesCount(garden->getBranchesTotal() + 1);
    }
}

void TREE_CLASS::fadeTree() {
    if(height == 0) return;
    BRANCH_CLASS *branch = branches;
    while(branch != NULL){
        branch->fadeBranch();
        branch = branch->getNext();
    }
    if(branches != NULL){
        if(branches->getHeight() == height){
            BRANCH_CLASS *toDelete = branches;
            toDelete->cutBranch(0);
            branches = branches->getNext();
            branchesCount--;
            if(garden != NULL)
                garden->setBranchesCount(garden->getBranchesTotal() - 1);
            delete toDelete;
        }
    }
    height--;
}

void TREE_CLASS::harvestTree(unsigned int weight) {
    BRANCH_CLASS *branch = branches;
    while(branch != NULL){
        branch->harvestBranch(weight);
        branch = branch->getNext();
    }
}

void TREE_CLASS::cutTree(unsigned int height) {
    if(height >= this->height) return;
    this->height = height;
    BRANCH_CLASS *branch = branches;
    while(branch != NULL){
        if(branch->getHeight() > height){
            BRANCH_CLASS *toDelete = branch;
            toDelete->cutBranch(0);
            branches = branches->getNext();
            branch = branch->getNext();
            branchesCount--;
            if(garden != NULL)
                garden->setBranchesCount(garden->getBranchesTotal() - 1);
            delete toDelete;
        } else
            branch = branch->getNext();
    }
}

void TREE_CLASS::cloneBranch(BRANCH_CLASS *obj) {
    if(branches == NULL) return;
    BRANCH_CLASS *branch = branches;
    BRANCH_CLASS *prevBranch = NULL;
    BRANCH_CLASS *prevClone = NULL;
    BRANCH_CLASS *toReplace = NULL;
    while(branch != NULL){
        if(branch->getLength() == 0) {
            prevClone = prevBranch;
            toReplace = branch;
        }
        prevBranch = branch;
        branch = branch->getNext();
    }
    if(toReplace != NULL){
        BRANCH_CLASS *clone = new BRANCH_CLASS(*obj);
        clone->setTree(this);
        clone->setHeight(toReplace->getHeight());
        clone->setNext(toReplace->getNext());
        fruitsWeight += clone->getWeightsTotal();
        fruitsCount += clone->getFruitsTotal();
        if(garden != NULL){
            garden->setWeight(garden->getWeightsTotal() + clone->getWeightsTotal());
            garden->setFruitCount(garden->getFruitsTotal() + clone->getFruitsTotal());
        }
        BRANCH_CLASS *toDelete;
        if(prevClone == NULL){
            toDelete = branches;
            branches = clone;
        } else{
            toDelete = toReplace;
            prevClone->setNext(clone);
        }
        delete toDelete;
    }
}

BRANCH_CLASS* TREE_CLASS::getBranchPointer(unsigned int height) {
    BRANCH_CLASS *branch = branches;
    while(branch != NULL){
        if(branch->getHeight() == height) return branch;
        branch = branch->getNext();
    }
    return branch;
}

//--- GARDEN_CLASS -----------------------------------------------------------------------------------------------------

GARDEN_CLASS::~GARDEN_CLASS() {
    while(trees != NULL){
        TREE_CLASS *toDelete = trees;
        trees = trees->getNext();
        delete toDelete;
    }
}

void GARDEN_CLASS::plantTree() {
    TREE_CLASS *newTree = new TREE_CLASS(0, this);
    treesCount++;
    if(trees == NULL) {
        trees = newTree;
        lastTree = trees;
    }
    else if(trees->getNumber() > 0){
        newTree->setNext(trees);
        trees->setPrev(newTree);
        trees = newTree;
    } else {
        TREE_CLASS *tree = trees;
        while(tree->getNext() != NULL){
            if(tree->getNext()->getNumber() > (tree->getNumber() + 1)){
                newTree->setNumber(tree->getNumber() + 1);
                newTree->setNext(tree->getNext());
                newTree->setPrev(tree);
                tree->getNext()->setPrev(newTree);
                tree->setNext(newTree);
                return;
            }
            tree = tree->getNext();
        }
        newTree->setNumber(tree->getNumber() + 1);
        tree->setNext(newTree);
        newTree->setPrev(tree);
        lastTree = newTree;
    }
}

void GARDEN_CLASS::extractTree(unsigned int num) {
    if(trees == NULL) return;

    TREE_CLASS *tree = trees;
    TREE_CLASS *prevTree = NULL;
    TREE_CLASS *prevDel = NULL;
    TREE_CLASS *toDel = NULL;
    while(tree != NULL) {
        if(tree->getNumber() == num){
            toDel = tree;
            prevDel = prevTree;
            break;
        }
        prevTree = tree;
        tree = tree->getNext();
    }

    if(toDel != NULL){
        fruitsWeight -= toDel->getWeightsTotal();
        fruitsCount -= toDel->getFruitsTotal();
        branchesCount -= toDel->getBranchesTotal();
        treesCount--;
        if(prevDel == NULL)
            trees = toDel->getNext();
        else
            prevDel->setNext(toDel->getNext());
        delete toDel;

    }
}
/*
void GARDEN_CLASS::extractTree(unsigned int num) {
    if(trees == NULL) return;

    TREE_CLASS *tree = lastTree;
    TREE_CLASS *nextTree = NULL;
    TREE_CLASS *nextDel = NULL;
    TREE_CLASS *toDel = NULL;
    while(tree != NULL) {
        if(tree->getNumber() == num){
            toDel = tree;
            nextDel = nextTree;
            break;
        }
        nextTree = tree;
        tree = tree->getPrev();
    }

    if(toDel != NULL){
        fruitsWeight -= toDel->getWeightsTotal();
        fruitsCount -= toDel->getFruitsTotal();
        branchesCount -= toDel->getBranchesTotal();
        treesCount--;
        if(nextDel == NULL) {  //usuwam ostatni
            lastTree = toDel->getPrev();
            lastTree->setNext(NULL);
        }
        else {
            nextDel->setPrev(toDel->getPrev());
            toDel->getPrev()->setNext(nextDel);
        }
        if(toDel == lastTree)
            lastTree = lastTree->getPrev();
        delete toDel;

    }
}*/

void GARDEN_CLASS::growthGarden() {
    TREE_CLASS *tree = trees;
    while(tree != NULL){
        tree->growthTree();
        tree = tree->getNext();
    }
}

void GARDEN_CLASS::fadeGarden() {
    TREE_CLASS *tree = trees;
    while(tree != NULL){
        tree->fadeTree();
        tree = tree->getNext();
    }
}

void GARDEN_CLASS::harvestGarden(unsigned int weight) {
    TREE_CLASS *tree = trees;
    while(tree != NULL){
        tree->harvestTree(weight);
        tree = tree->getNext();
    }
}

TREE_CLASS* GARDEN_CLASS::getTreePointer(unsigned int num) {
    TREE_CLASS* tree = trees;
    while(tree != NULL){
        if(tree->getNumber() == num) return tree;
        tree = tree->getNext();
    }
    return tree;
}

void GARDEN_CLASS::cloneTree(unsigned int num) {
    TREE_CLASS *toClone = getTreePointer(num);
    if(toClone == NULL) return;

    TREE_CLASS *newTree = new TREE_CLASS(*toClone);
    newTree->setGarden(this);
    treesCount++;
    branchesCount += newTree->getBranchesTotal();
    fruitsCount += newTree->getFruitsTotal();
    fruitsWeight += newTree->getWeightsTotal();
    if(trees == NULL){
        trees = newTree;
        lastTree = trees;
    } else if(trees->getNumber() > 0){

        trees->setPrev(newTree);
        newTree->setNext(trees);
        trees = newTree;
    } else {
        TREE_CLASS *tree = trees;
        while(tree->getNext() != NULL){
            if(tree->getNext()->getNumber() > (tree->getNumber() + 1)){
                newTree->setNumber(tree->getNumber() + 1);

                newTree->setNext(tree->getNext());
                newTree->setPrev(tree);
                tree->getNext()->setPrev(newTree);
                tree->setNext(newTree);
                return;
            }
            tree = tree->getNext();
        }
        newTree->setNumber(tree->getNumber() + 1);

        newTree->setPrev(tree);
        tree->setNext(newTree);
        lastTree = newTree;
    }
}

void AFR(int *, unsigned short int **, int ***);
void ALR(int *, unsigned short int **, int ***);
void AFC(int *, unsigned short int **, int ***);
void ALC(int *, unsigned short int **, int ***);
void IBR(int *, unsigned short int **, int ***);
void IAR(int *, unsigned short int **, int ***);
void IBC(int *, unsigned short int **, int ***);
void IAC(int *, unsigned short int **, int ***);
void SWR(int *, unsigned short int **, int ***);
void SWC(int *, unsigned short int **, int ***);
void DFR(int *, unsigned short int **, int ***);
void DLR(int *, unsigned short int **, int ***);
void DFC(int *, unsigned short int **, int ***);
void DLC(int *, unsigned short int **, int ***);
void RMR(int *, unsigned short int **, int ***);
void RMC(int *, unsigned short int **, int ***);
void RMB(int *, unsigned short int **, int ***);
void ISB(int *, unsigned short int **, int ***);
void WRF(int *, unsigned short int **, int ***);
void RDF(int *, unsigned short int **, int ***);
void PRT(int *, unsigned short int **, int ***);
void END(int *, unsigned short int **, int ***);

void callCommand(char*, int*, unsigned short int **, int ***);

void insertRow(int* , unsigned short int **, int ***, int, unsigned short int);
void insertValsInRow(unsigned short int **, int ***, int, int, int);
void removeRow(int *, unsigned short int **, int ***, int);
void removeValsInRow(unsigned short int **, int ***, int, int, int);

void insertRow(int* rowsNum, unsigned short int **columnsNum, int ***array, int idx, unsigned short int size){
    *array = (int**) realloc(*array, (*rowsNum + 1) * sizeof(int*));
    *columnsNum = (unsigned short int*) realloc(*columnsNum, (*rowsNum + 1) * sizeof(unsigned short int));
    int i;
    for(i = *rowsNum ; i > idx ; i--){
        *(*array + i) = *(*array + i - 1);
        *(*columnsNum + i) = *(*columnsNum + i - 1);
    }
    *rowsNum += 1;
    *(*columnsNum + idx) = size;
    *(*array + idx) = (int*) malloc(size * sizeof(int));
    for(i = 0 ; i < size ; i++){
        scanf("%d", (*(*array + idx) + i));
    }
}

void insertValsInRow(unsigned short int **columnsNum, int ***array, int rowIdx, int columnIdx, int n){
    *(*array + rowIdx) = (int*) realloc(*(*array + rowIdx), (*(*columnsNum + rowIdx) + n) * sizeof(int));
    int i;
    for(i = (*(*columnsNum + rowIdx) + n - 1) ; i >= columnIdx + n ; i--){
        *(*(*array + rowIdx) + i) = *(*(*array + rowIdx) + i - n);
    }
    *(*columnsNum + rowIdx) += n;
    for(i = columnIdx ; i < columnIdx + n ; i++) {
        scanf("%d", (*(*array + rowIdx) + i));
    }
}

void removeRow(int* rowsNum, unsigned short int **columnsNum, int ***array, int rowIdx){
    free(*(*array + rowIdx));
    * rowsNum -= 1;
    int i;
    for(i = rowIdx ; i < *rowsNum ; i++){
        *(*array + i) = *(*array + i + 1);
        *(*columnsNum + i) = *(*columnsNum + i + 1);
    }
    *(*array + i) = NULL;
    *array = (int**) realloc(*array, *rowsNum * sizeof(int*));
    *columnsNum = (unsigned short int*) realloc(*columnsNum, *rowsNum * sizeof(unsigned short int));
}

void removeValsInRow(unsigned short int **columnsNum, int ***array, int rowIdx, int columnIdx, int n){
    int i;
    for(i = columnIdx + n ; i < *(*columnsNum + rowIdx) ; i++){
        *(*(*array + rowIdx) + i - n) = *(*(*array + rowIdx) + i);
    }
    if(columnIdx + n >= *(*columnsNum + rowIdx)){
        *(*columnsNum + rowIdx) = columnIdx;
    }
    else{
        *(*columnsNum + rowIdx) -= n;
    }
    *(*array + rowIdx) = (int*) realloc((*(*array + rowIdx)), *(*columnsNum + rowIdx) * sizeof(int));
}

void AFR(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int w;
    scanf("%d", &w);
    if(w <= 0) return;
    insertRow(rowsNum, columnsNum, array, 0, w);
}

void ALR(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int w;
    scanf("%d", &w);
    if(w <= 0) return;
    insertRow(rowsNum, columnsNum, array, *rowsNum, w);
}

void AFC(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int h;
    scanf("%d", &h);
    if(h <= 0) return;
    int i;
    for(i = 0 ; i < *rowsNum && i < h; i++){
        insertValsInRow(columnsNum, array, i, 0, 1);
    }
    for(i ; i < h ; i++){
        insertRow(rowsNum, columnsNum, array, *rowsNum, 1);
    }
}

void ALC(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int h;
    scanf("%d", &h);
    if(h <= 0) return;
    int i;
    for(i = 0 ; i < *rowsNum && i < h; i++){
        insertValsInRow(columnsNum, array, i, *(*columnsNum + i), 1);
    }
    for(i ; i < h ; i++){
        insertRow(rowsNum, columnsNum, array, *rowsNum, 1);
    }
}

void IBR(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int r;
    scanf("%d", &r);
    if(r >= *rowsNum || r < 0) return;
    int w;
    scanf("%d", &w);
    if(w <= 0) return;
    insertRow(rowsNum, columnsNum, array, r, w);
}

void IAR(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int r;
    scanf("%d", &r);
    if(r >= *rowsNum || r < 0) return;
    int w;
    scanf("%d", &w);
    insertRow(rowsNum, columnsNum, array, r + 1, w);
}


void IBC(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int c, h;
    scanf("%d %d", &c, &h);
    if(h <= 0) return;
    int i;
    for(i = 0 ; i < *rowsNum && i < h; i++){///
        int idx = c;
        if(idx > *(*columnsNum + i)){
            idx = *(*columnsNum + i);
        }
        insertValsInRow(columnsNum, array, i, idx, 1);
    }
    for(i ; i < h ; i++){
        insertRow(rowsNum, columnsNum, array, *rowsNum, 1);
    }
}

void IAC(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int c, h;
    scanf("%d %d", &c, &h);
    if(h <= 0) return;
    c += 1;
    int i;
    for(i = 0 ; i < *rowsNum && i < h; i++){
        int idx = c;
        if(idx > *(*columnsNum + i)){
            idx = *(*columnsNum + i);
        }
        insertValsInRow(columnsNum, array, i, idx, 1);
    }
    for(i ; i < h ; i++){
        insertRow(rowsNum, columnsNum, array, *rowsNum, 1);
    }
}

void ISB(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int r, h, c, w;
    scanf("%d %d %d %d", &r, &c, &h, &w);
    if(r > *rowsNum) {
        r = *rowsNum;
    }
    h += r;
    int i;
    for(i = r ; i < h && i < *rowsNum ; i++){
        if(c > *(*columnsNum + i)) {
            insertValsInRow(columnsNum, array, i, *(*columnsNum + i), w);
        } else{
            insertValsInRow(columnsNum, array, i, c, w);
        }
    }
    for(i ; i < h ; i++){
        insertRow(rowsNum, columnsNum, array, *rowsNum, w);
    }
}

void SWR(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int r, s;
    scanf("%d %d", &r, &s);
    if(*rowsNum == 0) return;
    if(r >= *rowsNum || s >= *rowsNum) return;
    if(r == s) return;
    int* tmpRow = *(*array + r);
    *(*array + r) = *(*array + s);
    *(*array + s) = tmpRow;
    int tmp = *(*columnsNum + r);
    *(*columnsNum + r) = *(*columnsNum + s);
    *(*columnsNum + s) = tmp;
}

void SWC(int *rowsNum, unsigned short int **columnsNum, int ***array){
    int c, d;
    scanf("%d %d", &c, &d);
    if(*rowsNum == 0) return;
    if(c == d) return;
    int i;
    for(i = 0 ; i < *rowsNum ; i++){
        if(c < *(*columnsNum + i) && d < *(*columnsNum + i)){
            int tmp = *(*(*array + i) + c);
            *(*(*array + i) + c) = *(*(*array + i) + d);
            *(*(*array + i) + d) = tmp;
        }
    }
}

void DFR(int *rowsNum, unsigned short int **columnsNum, int ***array){
    if(*rowsNum == 0) return;
    removeRow(rowsNum, columnsNum, array, 0);
}

void DLR(int *rowsNum, unsigned short int **columnsNum, int ***array){
    if(*rowsNum == 0) return;
    removeRow(rowsNum, columnsNum, array, *rowsNum - 1);
}

void DFC(int *rowsNum, unsigned short int **columnsNum, int ***array){
    if(*rowsNum == 0) return;
    int i;
    for(i = 0 ; i < *rowsNum ; i++){
        if(*(*columnsNum + i) == 1){
            removeRow(rowsNum, columnsNum, array, i);
            i--;
        } else{
            removeValsInRow(columnsNum, array, i, 0, 1);
        }
    }
}

void DLC(int* rowsNum, unsigned short int **columnsNum, int ***array){
    if(*rowsNum == 0) return;
    int i;
    for(i = 0 ; i < *rowsNum ; i++){
        if(*(*columnsNum + i) == 1) {
            removeRow(rowsNum, columnsNum, array, i);
            i--;
        } else{
            removeValsInRow(columnsNum, array, i, *(*columnsNum + i) - 1, 1);
        }
    }
}

void RMR(int* rowsNum, unsigned short int **columnsNum, int ***array){
    int r;
    scanf("%d", &r);
    if(*rowsNum == 0) return;
    if(r >= *rowsNum) return;
    removeRow(rowsNum, columnsNum, array, r);
}

void RMC(int *rowsNum, unsigned short int **columnsNum, int ***array){
    int c;
    scanf("%d", &c);
    if(*rowsNum == 0) return;
    if(c == 0){
        DFC(rowsNum, columnsNum, array);
        return;
    }
    int i;
    for(i = 0 ; i < *rowsNum ; i++){
        if(c < *(*columnsNum + i)){
            removeValsInRow(columnsNum, array, i, c, 1);
        }
    }
}

void RMB(int *rowsNum, unsigned short int **columnsNum, int ***array){
    int r, h, c, w;
    scanf("%d %d %d %d", &r, &h, &c, &w);
    if(*rowsNum == 0) return;
    if(r >= *rowsNum) return;
    int i;
    for(i = r ; i < r + h && i < *rowsNum ; i++){
        if((c == 0) && (*(*columnsNum + i) <= w)){
            removeRow(rowsNum, columnsNum, array, i);
            i--;
            h--;
        } else{
            if(c < *(*columnsNum + i))
                removeValsInRow(columnsNum, array, i, c, w);
        }
    }
}

void PRT(int *rowsNum, unsigned short int **columnsNum, int ***array) {
    int i;
    int j;
    printf("%d\n", *rowsNum);
    for (i = 0; i < *rowsNum; i++) {
        printf("%d", *(*columnsNum + i));
        for (j = 0; j < *(*columnsNum + i); j++) {
            printf(" %d", *(*(*array + i) + j));
        }
        printf("\n");
    }
}

void END(int* rowsNum, unsigned short int **columnsNum, int ***array){
    if(*rowsNum > 0) {
        int i;
        for (i = 0; i < *rowsNum; i++) {
            free(*(*array + i));
        }
        free(*array);
        *array = NULL;
        free(*columnsNum);
        *columnsNum = NULL;
        *rowsNum = 0;
    }
    // jezeli rowsNum == 0 to pamiec w array i columnsNum juz byla zwolniona
}

unsigned short int endianUnsgnShortConversion(unsigned short int val){
    return ((val >> 8) | (val << 8));
}

int* endianIntArrConversion(int** valsArr, int size){
    int* arr = (int*)malloc(size * sizeof(int));
    int i;
    for(i = 0 ; i < size ; i++){
        *(arr + i) = ((*(*valsArr + i) >> 24) & 0x000000ff) | ((*(*valsArr + i) << 8) & 0x00ff0000) |
                     ((*(*valsArr + i) >> 8) & 0x0000ff00) | ((*(*valsArr + i) << 24) & 0xff000000);
    }
    return arr;
}

void WRF(int* rowsNum, unsigned short int **columnsNum, int ***array){
    char*  filename = (char*) malloc(100 * sizeof(char));
    scanf("%s", filename);
    FILE *myFile;
    myFile = fopen(filename, "wb");
    if(myFile != NULL) {
        fprintf(myFile, "%d\n", *rowsNum);
        int* row;
        unsigned short int columns;
        int i;
        for(i = 0 ; i < *rowsNum ; i++){
            columns = endianUnsgnShortConversion(*(*columnsNum + i));
            row = endianIntArrConversion((*array + i), *(*columnsNum + i));
            fwrite(&columns, sizeof(unsigned short int), 1, myFile);
            fwrite(row, sizeof(int), *(*columnsNum + i), myFile);
            free(row);
            row = NULL;
        }
        if (fclose(myFile) != 0){
            printf("Problem with closing the file %s\n", filename);
        }
    } else{
        printf("Problem with opening the file %s\n", filename);
    }
    free(filename);
}

void RDF(int* rowsNum, unsigned short int **columnsNum, int ***array){
    if(*rowsNum != 0)
        END(rowsNum, columnsNum, array);
    char*  filename = (char*) malloc(100 * sizeof(char));
    scanf("%s", filename);
    FILE *myFile;
    myFile = fopen(filename, "rb");
    if(myFile != NULL) {
        fscanf(myFile, "%d\n", rowsNum);
        *columnsNum = (unsigned short int*) malloc(*rowsNum * sizeof(unsigned short int));
        *array = (int**) malloc((*rowsNum) * sizeof(int*));
        int* row;
        unsigned short int columns;
        int i;
        for(i = 0 ; i < *rowsNum ; i++){
            fread(&columns,sizeof(unsigned short int), 1, myFile);
            *(*columnsNum + i) = endianUnsgnShortConversion(columns);
            row = (int*) malloc((*(*columnsNum + i)) * sizeof(int));
            fread(row, sizeof(int), *(*columnsNum + i), myFile);
            *(*array + i) = endianIntArrConversion(&row, *(*columnsNum + i));
            free(row);
            row = NULL;
        }
        if (fclose(myFile) != 0){
            printf("Problem with closing the file %s\n", filename);
        }
    } else{
        printf("Problem with opening the file %s\n", filename);
    }
    free(filename);
}

void callCommand(char* command, int* numOfRows, unsigned short int **sizeOfRows, int ***array2D) {
    if(strcmp(command, "\n") == 0) {  // fgets traktuje '\n' jako poprawny znak i go pobiera, wykluczam ten przypadek
        return;
    } else if(strcmp(command, "AFR") == 0) {
        AFR(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "ALR") == 0) {
        ALR(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "AFC") == 0) {
        AFC(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "ALC") == 0) {
        ALC(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "IBR") == 0) {
        IBR(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "IAR") == 0) {
        IAR(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "IBC") == 0) {
        IBC(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "IAC") == 0) {
        IAC(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "SWR") == 0) {
        SWR(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "SWC") == 0) {
        SWC(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "DFR") == 0) {
        DFR(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "DLR") == 0) {
        DLR(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "DFC") == 0) {
        DFC(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "DLC") == 0) {
        DLC(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "RMR") == 0) {
        RMR(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "RMC") == 0) {
        RMC(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "RMB") == 0) {
        RMB(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "ISB") == 0) {
        ISB(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "WRF") == 0) {
        WRF(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "RDF") == 0) {
        RDF(numOfRows, sizeOfRows, array2D);
    } else if(strcmp(command, "PRT") == 0) {
        PRT(numOfRows, sizeOfRows, array2D);
    }
}

int main() {

    matrix A = {{2, 1, 2}, {2, -1, 3}, {4, 6, 4}};
    vector b = {-3, 2, 1};
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "Rozwiazanie rownania Ax = b przy uzyciu eliminacji gaussa, "
                 "wykorzystując trzy metody wyboru elementu glownego:" << std::endl;
    std::cout << "---" << std::endl;
    std::cout << "A = " << std::endl;
    printMat(A);
    std::cout << "---" << std::endl;
    std::cout << "b = " << std::endl;
    printVec(b);
    std::cout << std::endl;
    std::cout << "---" << std::endl;

    std::cout << "wybor czesciowy: " << std::endl;
    vector x = Gauss(A, b, maxPartialChoice);
    printVec(x);
    std::cout << std::endl;
    std::cout << "---" << std::endl;
    std::cout << "wybor pelny: " << std::endl;
    x = Gauss(A, b, maxFullChoice);
    printVec(x);
    std::cout << std::endl;
    std::cout << "---" << std::endl;
    std::cout << "wybor skalowany: " << std::endl;
    x = Gauss(A, b, maxScaledPartialChoice);
    printVec(x);
    std::cout << std::endl;
    std::cout << "---" << std::endl;
    std::cout << "rozwiazanie dokladne:" << std::endl;
    std::cout << "-10.85 1.75 8.5" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    return 0;

    int **array2D = NULL;  // tablica 2D
    int numberOfRows = 0;  // przechowuje ilosc wierszy tablicy
    unsigned short int *sizesOfRows = NULL;  // przechowuje ilosc kolumn dla kazdego wiersza

    char* command;
    command = (char*) malloc(4 * sizeof(char));
    fgets(command, 4, stdin);  // pobiera pierwsze trzy znaki z stdin (komenda)

    while(strcmp(command, "END") != 0){
        callCommand(command, &numberOfRows, &sizesOfRows, &array2D);
        fgets(command, 4, stdin);
    }

    free(command);
    END(&numberOfRows, &sizesOfRows, &array2D);

//    return 0;

    Node* list = NULL;

    char command2;
    std::cin >> command2;
    while(command2 != 'F') {
        callCommand(command2, &list);
        std::cin >> command2;
    }
    clearList(&list);
    return 0;
}