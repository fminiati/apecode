#ifndef FILEPARSER_H
#define FILEPARSER_H

#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace fm
{
    class FileParser
    {
    public:
        // destructor
        ~FileParser() {}

        // constructor
        FileParser(const std::string a_file_name)
        {
            // open and parse file
            if (std::ifstream file(a_file_name); file)
            {
                bool append_entry_to_prev_key{false};

                std::string line, key;
                while (getline(file, line))
                {
                    // valid substr: ignorre initial blanks and commented out text
                    const auto beg = line.find_first_not_of(' ');
                    const auto end = line.find_first_of('#');
                    if (end > beg)
                    {
                        // stop short of backslash too, indicating input continues onto next line
                        const auto last{line.find_last_not_of("\\ #", end)};
                        // define string to work with
                        const std::string s{line.substr(beg, last - beg + 1)};

                        if (append_entry_to_prev_key)
                        {
                            // add to previous key
                            m_entries[key].append(" " + s);
                        }
                        else if (const auto pos = s.find('='); pos != std::string::npos)
                        {
                            // store key and entry separated by '=' sign
                            key = s.substr(0, s.find_last_not_of(' ', pos - 1) + 1);

                            m_entries[key] = s.substr(pos + 1);
                        }

                        // is there a backslash before commented out text ?
                        append_entry_to_prev_key = false;
                        auto i = 1 + last;
                        for (; i < line.size() && line[i] == ' '; ++i)
                        {}
                        if (i < line.size() && line[i] == '\\')
                            append_entry_to_prev_key = true;
                    }
                    else
                    {
                        if (append_entry_to_prev_key)
                        {
                            std::cerr << "\nexpecting input continuation from previous line..."
                                         " perhaps forgot open \\ in input file!\n\n";
                            // close file
                            file.close();
                            std::abort();
                        }
                    }
                }
                // close file
                file.close();
            }
            else
            {
                std::cerr << "FileParser::define: Error while opening file!\n";
                std::abort();
            }
        }

        // get item value of entry a_name
        template <class T>
        void get_item(T &a_val, const std::string a_name)
        {
            // search object a_name
            if (const auto it = m_entries.find(a_name); it != m_entries.end())
            {
                try
                {
                    std::istringstream line(it->second);
                    line >> a_val;
                }
                catch (std::ios_base::failure &)
                {
                    std::cout << "FileParser::get_item: error: entry was : "
                              << it->second << "\n";
                }
            }
            else
            {
                std::cerr << "FileParser::get_item: item " << a_name << " not found!\n";
                std::abort();
            }
        }

        // fill a_c container with (a_c.size) items from entry a_name
        template <class T>
        void get_items(T &a_c, const std::string a_name)
        {
            assert(a_c.size() > 0); //"FileParser::get_items: size of "+a_name+" must be >0");

            // search object a_name
            if (const auto it = m_entries.find(a_name.c_str()); it != m_entries.end())
            {
                try
                {
                    std::istringstream line(it->second);
                    for (auto &c : a_c)
                        line >> c;
                }
                catch (std::ios_base::failure &)
                {
                    std::cout << "FileParser::get_item: error: entry was : "
                              << it->second << "\n";
                }
            }
            else
            {
                std::cerr << "FileParser::get_item: item " << a_name << " not found!\n";
                std::abort();
            }
        }

        // flush entries
        void flush() const
        {
            for (auto [key, entry] : m_entries)
            {
                std::cout << "\"" << key << "\" -- \"" << entry << "\"" << '\n';
            }
        }

    protected:
        std::unordered_map<std::string, std::string> m_entries;
    };
}; // namespace fm

#endif
