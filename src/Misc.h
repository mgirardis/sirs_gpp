#ifndef _MISCDEF_H
#define _MISCDEF_H

#define _CRT_SECURE_NO_WARNINGS
#include <sys/stat.h>
#include <stdio.h>
#include <cmath>
#include <memory>
#include <time.h>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>
#include <string.h>
#include <string>
#include <vector>

namespace Misc
{
    struct AvgVar_t
    {
        double Average = 0.0;
        double Variance = 0.0;
    };

    template <typename EnumType>
    inline EnumType StrToEnum(std::string s, const char* ec[])
    {
        for (int t = (int)EnumType::begin + 1; t < (int)EnumType::end; t++)
        {
            if (s.compare(static_cast<std::string>(ec[t])) == 0)
            {
                return (EnumType)t;
            }
        }
        throw std::invalid_argument("Unknown value for enum");
    };

    template <typename EnumType>
    inline std::string EnumToStr(EnumType ev, const char* ec[])
    {
        for (int t = (int)EnumType::begin + 1; t < (int)EnumType::end; t++)
        {
            if (ev == (EnumType)t)
            {
                return static_cast<std::string>(ec[t]);
            }
        }
        throw std::invalid_argument("Unknown value for enum");
    };

    template <typename EnumType>
    inline std::string EnumStrList(const char* ec[])
    {
        std::string s = "";
        for (int t = (int)EnumType::begin + 1; t < (int)EnumType::end; t++)
        {
            s += EnumToStr((EnumType)t, ec) + "; ";
        }
        return s;
    };

    static inline AvgVar_t CalcAverageVariance(std::vector<double>& v)
    {
        AvgVar_t a = {};
        int i = 0, n = static_cast<int>(v.size());
        double mean = 0.0;
        double var = 0.0;
        while (i < n)
        {
            mean += v[i];
            var += v[i] * v[i];
            i++;
        }
        mean /= static_cast<double>(n);
        var = var / static_cast<double>(n) - mean * mean; //= sqrt(std / static_cast<double>(n) - mean * mean);
        a.Average = mean;
        a.Variance = (var < 0.0 ? 0.0 : var);
        return a;
    }

    static inline std::string GetSimulationTime(clock_t t0, clock_t t1)
    {
        double totalElapsedSeconds = double(t1 - t0) / CLOCKS_PER_SEC;
        int elapsedHours = static_cast<int>(floor(totalElapsedSeconds / 3600.0));
        int elapsedMinutes = static_cast<int>(floor(fmod(totalElapsedSeconds / 60.0, 60.0)));
        int elapsedSeconds = static_cast<int>(floor(fmod(fmod(totalElapsedSeconds, 60.0), 60.0)));
        int elapsedMilisec = static_cast<int>((totalElapsedSeconds - floor(totalElapsedSeconds)) * 1.0E3);

        char buffer[100];
#ifdef _WIN32
        sprintf_s(buffer, "Total simulation time (H:m:s.ms) = %02d:%02d:%02d.%03d", elapsedHours, elapsedMinutes, elapsedSeconds, elapsedMilisec);
#else
        sprintf(buffer, "Total simulation time (H:m:s.ms) = %02d:%02d:%02d.%03d", elapsedHours, elapsedMinutes, elapsedSeconds, elapsedMilisec);
#endif
        std::string res = std::string(buffer);
        return res;
    }

    static inline std::string RunProcess(std::string cmd){
#ifdef _WIN32
        FILE* pipe = _popen(cmd.c_str(), "r");
#else
        FILE* pipe = popen(cmd.c_str(), "r");
#endif
        if (!pipe) return "ERROR";
        char buffer[262144];
        std::string data;
        std::string result;
        int dist = 0;
        int size;
        //TIME_START
        while (!feof(pipe)) {
            size = (int)fread(buffer, 1, 262144, pipe); //cout<<buffer<<" size="<<size<<endl;
            data.resize(data.size() + size);
            memcpy(&data[dist], buffer, size);
            dist += size;
        }
        //TIME_PRINT_
#ifdef _WIN32
        _pclose(pipe);
#else
        pclose(pipe);
#endif
        return data;
    }

    static inline bool IsValidArgValue(std::string& str)
    {
        return str.find_first_not_of("+-.0123456789eE") == std::string::npos;
    }

    static inline std::string& StringLTrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    static inline bool StringContains(std::string& s1, std::string& s2)
    {
        if (s1.find(s2) != std::string::npos)
            return true;
        return false;
    }

    // trim from end
    static inline std::string& StringRTrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    // trim from both ends
    static inline std::string& StringTrim(std::string &s) {
        return Misc::StringLTrim(Misc::StringRTrim(s));
    }

    static inline std::string StringTrimAndReduce(std::string& input)
    {
        std::string output;
        std::unique_copy(input.begin(), input.end(), std::back_insert_iterator<std::string>(output),
            [](char a, char b){ return std::isspace(a) && std::isspace(b); });
        return Misc::StringTrim(output);
    }

    static inline std::string StringToLower(std::string str)
    {
        int length = str.size();
        int i = 0;
        while (i < length)
        {
            str[i] = tolower(str[i]);
            i++;
        }
        return str;
    }

    static inline std::vector<std::string> StringSplit(std::string str, std::string delimiter)
    {
        std::vector<std::string> arr;

        if (str.find_first_of(delimiter.c_str()) != -1)
        {
            int strleng = str.length();
            int delleng = delimiter.length();
            if (delleng == 0)
            {
                arr.push_back(str);//no change
                return arr;
            }

            int i = 0;
            int k = 0;
            while (i < strleng)
            {
                int j = 0;
                while (i + j < strleng && j < delleng && str[i + j] == delimiter[j])
                    j++;
                if (j == delleng)//found delimiter
                {
                    arr.push_back(str.substr(k, i - k));
                    i += (int)delleng;
                    k = i;
                }
                else
                {
                    i++;
                }
            }
            arr.push_back(str.substr(k, i - k));
            return arr;
        }
        else
        {
            arr.push_back(str);//no change
            return arr;
        }
    }

    static inline std::string StringJoin(std::vector<std::string> pieces, std::string glue)
    {
        std::string a;
        int leng = pieces.size();
        int i = 0;
        while (i < leng)
        {
            a += pieces[i];
            if (i < (leng - 1))
            {
                a += glue;
            }
            i++;
        }
        return a;
    }

    static inline bool FileExists(std::string fileName)
    {
        struct stat buffer;
        return (stat(fileName.c_str(), &buffer) == 0);
    }

    static inline std::string CheckAndGetFileName(std::string fileName, std::string fileExt)
    {
        std::string suffix = fileExt;
        char count[200];
        int i = 0;
        while (Misc::FileExists(fileName + suffix))
        {
#ifdef _WIN32
            sprintf_s(count, "_%d", ++i);
#else
            sprintf(count, "_%d", ++i);
#endif
            suffix = (std::string)count + fileExt;
        }
        return fileName + suffix;
    }

    static inline std::string CheckAndGetFileName(std::string fileName)
    {
        std::vector<std::string> v = StringSplit(fileName, ".");
        std::string str;
        std::string ext;
        if (v.size() > 1)
        {
            ext = "." + v.back();
            v.pop_back();
            str = StringJoin(v, ".");
        }
        else
        {
            str = fileName;
            ext = "";
        }
        return Misc::CheckAndGetFileName(str, ext);
    }

    static inline std::string GetFileNameWithoutExtension(std::string fileName)
    {
        std::vector<std::string> filenameParts;
        filenameParts = StringSplit(fileName, ".");
        if (filenameParts.size() > 1)
        {
            filenameParts.pop_back();
        }
        fileName = StringJoin(filenameParts, ".");
        return fileName;
    }
}

#endif
