'use strict';
const canvas = document.getElementById('gameCanvas');
const ctx = canvas.getContext('2d');
const grid = 20;
let count = 0;
let snake = { x: 160, y: 160, cells: [], maxCells: 4 };
let apple = { x: 320, y: 320 };
let dx = grid;
let dy = 0;

function getRandomInt(min, max) {
  return Math.floor(Math.random() * (max - min)) + min;
}

document.addEventListener('keydown', (e) => {
  switch (e.code) {
    case 'ArrowLeft':
      if (dx === 0) { dx = -grid; dy = 0; }
      break;
    case 'ArrowUp':
      if (dy === 0) { dx = 0; dy = -grid; }
      break;
    case 'ArrowRight':
      if (dx === 0) { dx = grid; dy = 0; }
      break;
    case 'ArrowDown':
      if (dy === 0) { dx = 0; dy = grid; }
      break;
  }
});

function loop() {
  requestAnimationFrame(loop);

  if (++count < 4) return;
  count = 0;

  ctx.clearRect(0, 0, canvas.width, canvas.height);

  snake.x += dx;
  snake.y += dy;

  if (snake.x < 0) {
    snake.x = canvas.width - grid;
  } else if (snake.x >= canvas.width) {
    snake.x = 0;
  }
  if (snake.y < 0) {
    snake.y = canvas.height - grid;
  } else if (snake.y >= canvas.height) {
    snake.y = 0;
  }

  snake.cells.unshift({ x: snake.x, y: snake.y });
  if (snake.cells.length > snake.maxCells) {
    snake.cells.pop();
  }

  ctx.fillStyle = 'red';
  ctx.fillRect(apple.x, apple.y, grid - 1, grid - 1);

  ctx.fillStyle = 'green';
  snake.cells.forEach((cell, index) => {
    ctx.fillRect(cell.x, cell.y, grid - 1, grid - 1);

    if (cell.x === apple.x && cell.y === apple.y) {
      snake.maxCells++;
      apple.x = getRandomInt(0, 20) * grid;
      apple.y = getRandomInt(0, 20) * grid;
    }

    for (let i = index + 1; i < snake.cells.length; i++) {
      if (cell.x === snake.cells[i].x && cell.y === snake.cells[i].y) {
        snake.x = 160;
        snake.y = 160;
        snake.cells = [];
        snake.maxCells = 4;
        dx = grid;
        dy = 0;
        apple.x = getRandomInt(0, 20) * grid;
        apple.y = getRandomInt(0, 20) * grid;
      }
    }
  });
}

requestAnimationFrame(loop);
