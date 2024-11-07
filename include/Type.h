#ifndef TYPEHOLDER_H
#define TYPEHOLDER_H
#include <cstring>

template<typename T>
class type {  // Изменено имя класса на 'type'
public:
    // Конструктор
    type(T value = 0.0) : value(value) {}

    // Оператор вывода
    friend std::ostream& operator<<(std::ostream& os, const type& t) {
        os << t.value;  // Вывести значение
        return os;
    }

    // Оператор ввода
    friend std::istream& operator>>(std::istream& is, type& t) {
        is >> t.value;  // Ввести значение
        return is;
    }

    // Переопределение арифметических операторов
    type operator+(const type& other) const {
        return type(value + other.value);
    }
    
    type operator-(const type& other) const {
        return type(value - other.value);
    }

    type operator*(const type& other) const {
        return type(value * other.value);
    }

    type operator/(const type& other) const {
        if (other.value != 0) {
            return type(value / other.value);
        }
        throw std::invalid_argument("Division by zero");
    }

    // Переопределение оператора присваивания
    type& operator=(const type& other) {
        if (this != &other) {
            value = other.value;
        }
        return *this;
    }

    // Переопределение оператора сравнения
    bool operator==(const type& other) const {
        return value == other.value;
    }
    T operator()(T val){
        if(std::strcmp(typeid(T).name(),"f")==0)return 1.19E-7;
        return 2.20E-16;
    }
private:
    T value;  // Хранит значение типа double
};

#endif