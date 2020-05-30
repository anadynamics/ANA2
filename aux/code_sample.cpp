#include <iterator>
#include <algorithm>
#include <utility>
#include <string_view>

template<typename Iterator, typename UnaryPredicate>
[[nodiscard]] constexpr auto count_if_to_last(Iterator begin, Iterator end, UnaryPredicate predicate){
    std::reverse_iterator rbegin{end};
    const std::reverse_iterator rend{begin};

    const auto last_item = std::find_if(rbegin, rend, predicate);
    const auto count = std::count_if(last_item, rend, predicate);

    return std::pair{count, last_item.base()};
}

template<typename Iterator, typename Value>
[[nodiscard]] constexpr auto count_to_last(Iterator begin, Iterator end, const Value &value)
{
    return count_if_to_last(begin, end, [&](const auto &input){ return input == value; });
}

int main(int argc, const char *argv[])
{
    constexpr std::string_view input{"hello\nworld"};
    constexpr auto result = count_to_last(input.begin(), input.end(), '\n');

    static_assert(result.first == 1);
    static_assert(result.second == std::next(input.begin(), 6));

    std::string_view runtime_input(argv[0]);
    return count_to_last(runtime_input.begin(), runtime_input.end(), '_').first;
}

