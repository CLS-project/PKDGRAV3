// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#pragma once

#include <typeindex>
#include <utility>
#include <any>
#include <boost/unordered/unordered_flat_map.hpp>

// The PKD, PST and LCL classes have a tendency to accumulate data that is used by multiple related services.
// To avoid bloating these classes with a large number of member variables that are only used by a few services,
// we provide a DataRegistry that allows services to store and retrieve data by type.
//
// Example usage:
//  struct MyServiceData {
//      int someValue;
//      std::vector<double> moreData;
//  };
//  // Storing data:
//  auto& data = registry.emplace<MyServiceData>();
//  data.someValue = 42;
//  data.moreData = {1.0, 2.0, 3.0};
//  // Retrieving data:
//  auto& data = registry.get<MyServiceData>();
//  std::cout << data.someValue << std::endl;
class DataRegistry {
  public:
    DataRegistry() = default;

    // Insert or replace
    template <typename T>
    auto set(T value) -> T& {
        auto& slot = data_[std::type_index(typeid(T))];
        slot = std::move(value);
        return std::any_cast<T&>(slot);
    }

    // Insert or replace data element (always constructs in-place)
    template <typename T, typename... Args>
    auto emplace(Args&&... args) -> T& {
        auto [it, inserted] = data_.try_emplace(std::type_index(typeid(T)), std::in_place_type<T>, std::forward<Args>(args)...);
        if (!inserted) it->second.template emplace<T>(std::forward<Args>(args)...);
        return std::any_cast<T&>(it->second);
    }
    // Retrieve or insert data element (construct in-place if missing)
    template <typename T, typename... Args>
    auto get_or_emplace(Args&&... args) -> T& {
        auto [it, inserted] = data_.try_emplace(std::type_index(typeid(T)), std::in_place_type<T>, std::forward<Args>(args)...);
        return std::any_cast<T&>(it->second);
    }

    // Access (throws if missing or wrong type)
    template <typename T>
    [[nodiscard]] auto get() -> T& {
        return std::any_cast<T&>(data_.at(std::type_index(typeid(T))));
    }
    template <typename T>
    [[nodiscard]] auto get() const -> const T& {
        return std::any_cast<const T&>(data_.at(std::type_index(typeid(T))));
    }

    // Query
    template <typename T>
    [[nodiscard]] auto contains() const noexcept -> bool {
        return data_.contains(std::type_index(typeid(T)));
    }

    // Remove
    template <typename T>
    void erase() noexcept {
        data_.erase(std::type_index(typeid(T)));
    }

    // Remove all
    void clear() noexcept { data_.clear(); }

  private:
    boost::unordered_flat_map<std::type_index, std::any> data_;
};

class DataRegistryOwner {
    DataRegistry data_;

  protected:
    // Retrieve data of type T. Throws if not present.
    template <typename T>
    auto data() -> T& {
        return data_.get<T>();
    }

    // Retrieve const data of type T. Throws if not present.
    template <typename T>
    auto data() const -> const T& {
        return data_.get<T>();
    }

    // Check if data of type T is present.
    template <typename T>
    [[nodiscard]] auto data_present() const noexcept -> bool {
        return data_.contains<T>();
    }

    // Remove data of type T.
    template <typename T>
    void data_erase() noexcept {
        data_.erase<T>();
    }

    // Insert or replace data of type T.
    template <typename T, typename... Args>
    auto data_emplace(Args&&... args) -> T& {
        return data_.emplace<T>(std::forward<Args>(args)...);
    }

    // Retrieve data of type T, constructing it in-place if not present.
    template <typename T, typename... Args>
    auto data_ensure(Args&&... args) -> T& {
        return data_.get_or_emplace<T>(std::forward<Args>(args)...);
    }
};
