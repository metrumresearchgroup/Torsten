-- Adds a unique link to each index entry in the HTML output
-- based on those special comments we have for functions

-- This is combined with the output of the gen_index.py script
-- to create a clickable index of all the functions in the documentation


index = quarto.utils.resolve_path("./functions-reference/functions_index.qmd")
indexText = io.open(index):read("*a")

function extractIndexEntry(elementText)
  if elementText:find("%; %-%-%>$") ~= nil then
    return "index-entry-" .. tostring(pandoc.sha1(elementText))
  end
  return nil
end

function escape(text)
  return text:gsub("\\", "\\\\"):gsub("*", "\\*")
end

if quarto.doc.is_format("html") then -- latex uses mkindex, not this
  return {
    RawBlock = function(el)
      -- this filter is in charge of producing the HTML anchors the index will link to
      if el.format == "html" then
        local indexEntry = extractIndexEntry(el.text)
        if indexEntry ~= nil then
          return pandoc.RawInline("html", "<a id=\"" .. indexEntry .. "\" class=\"index-entry\"></a> " .. el.text)
        end
      end
      return nil -- no change
    end,
    Strong = function(el2)
      return pandoc.walk_inline(el2, {
        Code = function(el3)
          -- this filter produces links to the index from the function's name
          -- Because we format these as **``functionName``** in the markdown
          -- it will always be a Strong element containing a Code element
          if el3.text ~= nil then
            -- only create a link if this appears in the index
            local escaped = escape(el3.text)
            if indexText:find(escaped, 1, true) ~= nil then
              return pandoc.Link(el3, "./functions_index.qmd#" .. escaped,
                "Jump to index entry", {class="unlink"})
            end
          end
          return nil
        end
      })
    end
  }

end
