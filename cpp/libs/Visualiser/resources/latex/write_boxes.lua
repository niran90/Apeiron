module(...,package.seeall)

-- Local constants
local show_boxes = false
local scale_dims = false
local scale      = 65536  -- 2^16sp = 1pt
local HLIST = node.id("hlist")
local VLIST = node.id("vlist")
local GLUE  = node.id("glue")
local KERN  = node.id("kern")
local GLYPH = node.id("glyph")


-- Helper functions
local function getgluespec(anode)
  if node.has_field(anode,"spec") then
    return anode.spec
  else
    return anode
  end
end

local function write_pos()
    local x, y = pdf.getpos()
    file = io.open("positions.txt", "a")
    if scale_dims then
      file:write(string.format("\n%.5f\t%.5f", x / scale, y / scale))
    else
      file:write('\n', x, '\t', y)
    end
    file:close()
end


-- The argument must be a box (hbox or vbox)
local function draw_elements(box)
  local parent = box
  local head   = box.list

  while head do
    if head.id == HLIST or head.id == VLIST then
      draw_elements(head)

      if show_boxes then
        local wd = head.width                  / scale
        local ht = (head.height + head.depth)  / scale
        local dp = head.depth                  / scale

        local pdfliteral = node.new("whatsit","pdf_literal")

        if head.id == HLIST then -- hbox
          pdfliteral.data = string.format("q 0.5 G 0.1 w 0 %g %g %g re s Q", -dp, wd, ht)
        else -- vbox
          pdfliteral.data = string.format("q 0.1 G 0.1 w 0 %g %g %g re s Q", 0, wd, -ht)
        end

        head.list = node.insert_before(head.list, head.list, pdfliteral)
      end

    elseif head.id == GLYPH then
      local c = string.utfcharacter(head.char)
      local wd = head.width  / scale
      local ht = head.height / scale
      local dp = head.depth  / scale

      file = io.open("attributes.txt", "a")
      if scale_dims then
        file:write('\n', c, string.format("\t%.5f\t%.5f\t%.5f", wd, ht, dp))
      else
        file:write('\n', c, '\t', head.width, '\t', head.height, '\t', head.depth)
      end
      file:close()

      local writepos = node.new("whatsit", "late_lua")
      writepos.data = write_pos
      -- _, head = node.insert_after(parent, head, writepos)
      parent.list, _ = node.insert_before(parent.list, head, writepos)

      if show_boxes then
        local pdfl = node.new("whatsit","pdf_literal")
        pdfl.data = string.format("q 0.7 G 0.1 w 0 %g %g %g re S f Q", -dp, -wd, dp + ht, -wd)
        _, head = node.insert_after(parent, head, pdfl)
      end

    elseif head.id == GLUE and show_boxes then
      local spec = getgluespec(head)
      local wd = spec.width -- the natural width of the glue
      local color = "0.5 G"
      
      if parent.glue_sign == 1 and parent.glue_order == spec.stretch_order then
        wd = wd + parent.glue_set * spec.stretch
        color = "0 0 1 RG"
      elseif parent.glue_sign == 2 and parent.glue_order == spec.shrink_order then
        wd = wd - parent.glue_set * spec.shrink
        color = "1 0 1 RG"
      end

      pdfliteral = node.new("whatsit","pdf_literal")

      if parent.id == HLIST then
        pdfliteral.data = string.format("q %s [0.2] 0 d  0.5 w 0 0  m %g 0 l s Q", color, wd / scale)
      else -- VLIST
        pdfliteral.data = string.format("q 0.1 G 0.1 w -0.5 0 m 0.5 0 l -0.5 %g m 0.5 %g l s [0.2] 0 d  0.5 w 0.25 0  m 0.25 %g l s Q",-wd / scale,-wd / scale,-wd / scale)
      end
      node.insert_before(parent.list,head,pdfliteral)

    else
      -- Any other node. Probably not interesting.
    end
    -- next node in our list. If the list is at the end, head becomes nil and
    -- the loop ends.
    head = head.next
  end

  return true
end


function shipout()
  -- Flush files
  io.open("positions.txt",  "w"):close()
  io.open("attributes.txt", "w"):close()

  -- StartWrite shipout box dimensions
  local shipout_box = tex.getbox("ShipoutBox")
  file = io.open("shipout_box.txt", "w")
  if scale_dims then
    file:write(string.format("%.5f\t%.5f\n", shipout_box.width / scale, shipout_box.height / scale))
  else
    file:write(shipout_box.width, "\t", shipout_box.height, "\n")
  end
  file:close()

  -- Recursively write and/or draw all glyph boxes
  draw_elements(shipout_box)
end